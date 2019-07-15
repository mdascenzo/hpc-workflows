# vi:syntax=python
import os, sys
import glob
from os import path
from snakemake.io import protected, expand, unpack
from snakemake.io import directory
from snakemake import shell, rules
from snakemake.utils import R
from snakemake.utils import update_config, available_cpu_count
from snakemake.exceptions import WorkflowError

# r2py reports all stdout using warnings
# ignore warnings to maintain cleaner stdout
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

# todo:
# 	- validate config input to ensure all files exist, generalize (e.g. all genome files, annotation files, etc.)
# 	- consider adding read length auto-detect
# 	- further generalize configuration section

# ---------------
# Configuration
# ---------------

# set default configuration file.
# also settable via commandline --configfile parameter which takes precedence over this assignment
if os.path.exists('config.yml'):
	configfile: 'config.yml'
if not config:
	raise WorkflowError('config file not found')

# load analysis options
opt_star = config['star']
opt_salmon = config['salmon']
# require that trim option be explicitly set in config, false otherwise
opt_trim = config['trim'] if 'trim' in config else False

SAMPLES = config['samples']
SAMPLE_IDS = list(SAMPLES.keys())

# set default options
opts = dict()
opts['salmon-quant-l'] = 'ISR'
opts['trimmomatic-adapters-fa'] = 'TruSeq3-PE-2.fa'
# override default options with config based options
if 'options' in config:
	opts.update(config['options'])

# trimmomatic.adapters_fa
# if a full path is not provided via opts use share directory relative to trimmomatic bin directory
if not os.path.isabs(opts['trimmomatic-adapters-fa']):
	opts['trimmomatic-adapters-fa'] = os.path.abspath(
		os.path.join(
			os.popen("echo `which trimmomatic`").read().rstrip(),
			'../../share/trimmomatic/adapters',
			opts['trimmomatic-adapters-fa']
		)
	)

# record versions of aligners used for index building purposes
star_version = os.popen("echo `star --version`").read().rstrip()
salmon_version = os.popen("echo `salmon --version`").read().rstrip().split(' ')[1]

# location of salmon index directory
salmon_index_location =\
	os.path.join(
		config['resources_dir'], 'transcriptomes', config['build'],
		config['tx_uid'], 'indexes/salmon', salmon_version
	) + '/'
# location of star index directory
star_index_location =\
	os.path.join(
		config['resources_dir'], 'genomes', config['build'],
		config['genome_uid'], 'indexes/star', star_version + '_sjo' + str(config['star_sj_db_overhang'])
	) + '/'


def set_log(name):
	return os.path.join(config['out'], name)

# ----------------
# workflow rules
# ----------------

rule all:
	input:
		# include fastqc
		expand(path.join(config['out'], 'fastqc/{sample}'), sample=SAMPLES),
		# include multiQC
		path.join(config['out'], "multiqc_report.html"),

		# option: include Salmon and tximport
		(
			path.join(config['out'], 'salmon', config['analysis_name'] + '_counts.txt'),
			path.join(config['out'], 'salmon', config['analysis_name'] + '_txi.rds'),
		) if opt_salmon else (),
		# option: include STAR and featureCounts
		(
			# star
			expand(
				path.join(config['out'], 'star/{sample}/Aligned.out.bam'),
				sample=SAMPLES
			),
			# featureCounts
			expand(
				path.join(config['out'], 'star/{sample}/feature_counts.txt'),
				sample=SAMPLES
			)

		) if opt_star else ()


def trim_reads_input(w):
	"""
	Evaluate global and local trim settings for processing by rule.trim_reads

	:param w: snakemake wildcard
	:return: input parameters
	"""

	input_parameters = {}

	if opt_trim or config['samples'][w.sample].get('trim'):

	 	input_parameters['read1'] = config['samples'][w.sample]['read1']
	 	input_parameters['read2'] = config['samples'][w.sample]['read2']

	return input_parameters


# rule.trim_reads
#
#	Description: 	Optionally trim adapter sequence from paired-end reads using Trimmomatic.
#
# 	Reference: 		Bolger, A.M., Lohse, M., Usadel, B., 2014. Trimmomatic: a flexible trimmer for Illumina sequence data.
# 					Bioinformatics 30, 2114–2120. https://doi.org/10.1093/bioinformatics/btu170
#
rule trim_reads:
	input:
		unpack(trim_reads_input)
	output:
		read1_paired = path.join(config['out'], 'trim_reads', '{sample}_R1_trimmed_paired.fq.gz'),
		read1_unpaired = path.join(config['out'], 'trim_reads', '{sample}_R1_trimmed_unpaired.fq.gz'),
		read2_paired =  path.join(config['out'], 'trim_reads', '{sample}_R2_trimmed_paired.fq.gz'),
		read2_unpaired = path.join(config['out'], 'trim_reads', '{sample}_R2_trimmed_unpaired.fq.gz')
	params:
		adapters = opts['trimmomatic-adapters-fa'],
		threads = available_cpu_count()-2
	shell:
		"""
		trimmomatic PE \
			-threads {params.threads} \
			{input.read1} \
			{input.read2} \
			{output.read1_paired} \
			{output.read1_unpaired} \
			{output.read2_paired} \
			{output.read2_unpaired} \
			ILLUMINACLIP:{params.adapters}:2:30:10:2:TRUE \
			MINLEN:25			
		"""


if opt_salmon:

	# rule.salmon_index
	#
	#	Description: 	Create index for Salmon quasi-mapping analysis
	#
	# 	Reference: 		Patro, R., Duggal, G., Love, M.I., Irizarry, R.A., Kingsford, C., 2017. Salmon provides fast and bias-aware
	# 					quantification of transcript expression.
	# 					Nature Methods 14, 417–419. https://doi.org/10.1038/nmeth.4197
	#
	#	Documentation:	https://salmon.readthedocs.io/en/latest/salmon.html
	#
	rule salmon_index:
		input:
			tx_fa = glob.glob(path.join(
				config['resources_dir'], 'transcriptomes', config['build'], config['tx_uid'], 'fa', '*.fa'
			))
		# index transcript fasta file if the directory does not already exist
		output:
			protected(
				directory(
					salmon_index_location
				)
			)
		log: 'salmon/salmon_index.log'

		threads: available_cpu_count()-2

		shell:
			"""
			salmon index --threads {threads} -t {input.tx_fa} -i {output} &> {log}
			"""

	def salmon_input(w):
		"""
		Provide rule.salmon_quant with either trimmed or un-trimmed reads given global and local trim settings

		:param w: snakemake wildcard
		:return: input parameters
		"""
		input = {}

		if opt_trim or config['samples'][w.sample].get('trim'):
			input['read1'] = rules.trim_reads.output.read1_paired
			input['read2'] = rules.trim_reads.output.read2_paired
		else:
			input['read1'] = config['samples'][w.sample]['read1']
			input['read2'] = config['samples'][w.sample]['read2']

		return input

	# rule.salmon_quant
	#
	# 	Reference:		Patro, R., Duggal, G., Love, M.I., Irizarry, R.A., Kingsford, C., 2017. Salmon provides fast and bias-aware
	# 					quantification of transcript expression.
	# 					Nature Methods 14, 417–419. https://doi.org/10.1038/nmeth.4197
	#
	#	Documentation:	https://salmon.readthedocs.io/en/latest/salmon.html
	#
	# 	Notes: 			Salmon libtype -  https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype
	#
	rule salmon_quant:
		input:
			unpack(salmon_input),
			index = salmon_index_location,
		output:
			quant = path.join(config['out'], 'salmon/{sample}/quant.sf')
		params:
			l = opts['salmon-quant-l']
		log: set_log('salmon/{sample}/stdout.log')

		threads: available_cpu_count()-2

		run:
			output_dir = os.path.dirname(output[0])
			cmd =\
			"""
			salmon quant \
				--threads {threads} \
				-i {input.index} \
				-l {params.l} \
				-1 {input.read1} \
				-2 {input.read2} \
				--validateMappings \
				--gcBias \
				--seqBias \
				-o {output_dir} &> {log}
			"""
			shell(cmd)



	# rule.aggrigate_tx
	#
	#	Description: 	Aggrigate transcript-level qunatification to gene-level.
	#
	# 	Reference:		Soneson, C., Love, M.I., Robinson, M.D., 2016. Differential analyses for RNA-seq:
	# 					transcript-level estimates improve gene-level inferences.
	# 					F1000Res 4, 1521–1521. https://doi.org/10.12688/f1000research.7563.2
	#
	#	Documentation:	https://bioconductor.org/packages/devel/bioc/manuals/tximport/man/tximport.pdf
	#	R-Vignette: 	https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
	#
	rule aggrigate_tx:
		input:
			tx2gene_fp = config['tx2gene_fp'],
			quant_files = expand(path.join(config['out'], 'salmon/{sample}/quant.sf'), sample=SAMPLES)
		output:
			counts = path.join(config['out'], 'salmon', config['analysis_name'] + '_counts.txt'),
			txi = path.join(config['out'], 'salmon', config['analysis_name'] + '_txi.rds')
		params:
			basename = config['analysis_name'],
			path = path.join(config['out'], 'salmon'),
			salmon_quant_analysis_dir = path.join(config['out'], 'salmon'),
			sample_ids = SAMPLE_IDS

		run:
			# noinspection RUnresolvedReference
			R("""
			source('./R/rnaseq_tools.R')
					
			# create R variables
			salmon_quant_analysis_dir <- "{params.salmon_quant_analysis_dir}" 
			sample_ids <- unlist(strsplit("{params.sample_ids}", " "))
			tx2gene_fp = "{input.tx2gene_fp}"
			basename = "{params.basename}"
			path = "{params.path}"
			
			# call rnaseq_tools::tximport.salmon()
			txi_salmon = tximport.salmon(salmon_quant_analysis_dir, sample_ids, tx2gene_fp)
			
			# save tx data to file
			write.txi(txi_salmon, basename=basename, path=path)
			
			""")

# STAR provides referenced based quantification:
# reference:
#
# optionally include STAR and downstream analysis components
#
if opt_star:


	# todo: check readlength, dynamically set sjdbOverhang
	# todo: allow mixed RL to be run in same analysis, currenly assumes the same for all samples
	# todo: create subdirectory based on GTF file, link or include GTF file with index

	# rule.salmon_index
	#
	#	Description: 	Create index for Salmon quasi-mapping analysis
	#
	# 	Reference: 		Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M.,
	# 					Gingeras, T.R., 2013. STAR: ultrafast universal RNA-seq aligner.
	# 					Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635
	#
	#	Documentation:	https://github.com/alexdobin/STAR/blob/2.6.1d/doc/STARmanual.pdf
	#
	#	Notes:			--sjdbOverhang: int>0: length of the donor/acceptor sequence on each side of the junctions
	#									Set to ReadLength-1 (e.g. 75-1 for 75bp reads)
	#
	rule star_index:
		input:
			fasta_files = path.join(
				config['resources_dir'], 'genomes', config['build'], config['genome_uid'], 'fa'
			),
			annotation_gtf = glob.glob(path.join(
				config['resources_dir'], 'genomes', config['build'], config['genome_uid'], 'annotation', '*.gtf'
			))[0]
		output:
			path = protected(
				directory(
					star_index_location
				)
			)
		params:
			sj_db_overhang = str(config['star_sj_db_overhang'])

		threads: available_cpu_count()-2

		run:
			command =\
			"star --runMode genomeGenerate --runThreadN {threads} --genomeFastaFiles {input.fasta_files}/*.fa --genomeDir {output.path} --outFileNamePrefix {output.path}"
			# add optional parameters
			if input.annotation_gtf is not None and config['star_sj_db_overhang'] is not None:
				command += ' --sjdbGTFfile ' + input.annotation_gtf
				command += ' --sjdbOverhang ' + params.sj_db_overhang

			shell(command)


	def star_align_input(w):
		"""
		Provide rule.star_align with either trimmed or un-trimmed reads given global and local trim settings

		:param w: snakemake wildcard
		:return: input parameters
		"""
		input = {}

		if opt_trim or config['samples'][w.sample].get('trim'):
			input['read1'] = rules.trim_reads.output.read1_paired
			input['read2'] = rules.trim_reads.output.read2_paired
		else:
			input['read1'] = config['samples'][w.sample]['read1']
			input['read2'] = config['samples'][w.sample]['read2']

		return input


	# rule.star_align
	#
	#	Description: 	Align paired-end reads
	#
	# 	Reference: 		Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M.,
	# 					Gingeras, T.R., 2013. STAR: ultrafast universal RNA-seq aligner.
	# 					Bioinformatics 29, 15–21. https://doi.org/10.1093/bioinformatics/bts635
	#
	#	Documentation:	https://github.com/alexdobin/STAR/blob/2.6.1d/doc/STARmanual.pdf
	#
	#	Notes: 			--genomeLoad: LoadAndKeep option appears to be osx incompatible, possibly containerize for improved performance
	#
	#					- Unsorted output can be directly input into featureCounts
	#					- Setting threads too high for this process may result in a ulimit error on osx
	rule star_align:
		input:
			unpack(star_align_input),
			star_genome_dir = rules.star_index.output.path
		output:
			bam = path.join(config['out'], 'star/{sample}/Aligned.out.bam'),
			count_file = path.join(config['out'], 'star/{sample}/ReadsPerGene.out.tab')
		params:
			quant_mode = 'TranscriptomeSAM GeneCounts',
			out_sam_type = 'BAM Unsorted'

		threads: 6

		shell:
			"""
			star --runThreadN {threads} \
				 --genomeDir {input.star_genome_dir} \
				 --readFilesIn <(gunzip -c {input.read1}) <(gunzip -c {input.read2}) \
				 --outFileNamePrefix  $(dirname {output.bam})/ \
				 --quantMode {params.quant_mode} \
				 --outSAMtype {params.out_sam_type}
			"""

	# rule.feature_counts
	#
	#	Description: 	Quantify alignments at gene-level
	#
	# 	Reference: 		Liao, Y., Smyth, G.K., Shi, W., 2014. featureCounts: an efficient general purpose program for assigning
	# 					sequence reads to genomic features.
	# 					Bioinformatics 30, 923–930. https://doi.org/10.1093/bioinformatics/btt656
	#
	#	Documentation:	https://github.com/alexdobin/STAR/blob/2.6.1d/doc/STARmanual.pdf
	#
	#	Notes:			-p include if read are paired-end
	#					-B only count read-pairs that have both ends aligned
	#					-s 0|1|2 : unstranded|stranded|reversely stranded
	#
	rule feature_counts:
		input:
			bam = rules.star_align.output.bam,
			annotation_gtf = glob.glob(path.join(
				config['resources_dir'], 'genomes', config['build'], config['genome_uid'], 'annotation', '*.gtf'
			))
		output:
			basename = path.join(config['out'], 'star/{sample}/feature_counts.txt'),
			summary = path.join(config['out'], 'star/{sample}/feature_counts.txt.summary')
		params:
			strandedness = 2

		threads: available_cpu_count()-2

		shell:
			"""
			featureCounts -T {threads} -p -B -s {params.strandedness} -a {input.annotation_gtf} -o {output.basename} {input.bam}
			"""


def fastqc_input(w):
	"""
	Provide rule.fastqc with either trimmed or un-trimmed reads given global and local trim settings

	:param w: snakemake wildcard
	:return: input parameters
	"""
	input = {}

	if opt_trim or config['samples'][w.sample].get('trim'):
		input['read1'] = rules.trim_reads.output.read1_paired
		input['read2'] = rules.trim_reads.output.read2_paired
	else:
		input['read1'] = config['samples'][w.sample]['read1']
		input['read2'] = config['samples'][w.sample]['read2']

	return input

# rule.fastqc
#
#	Description: 	Provides general QC metrics
#
# 	Reference: 		Andrews S. (2010). FastQC: a quality control tool for high throughput sequence data.
# 					http://www.bioinformatics.babraham.ac.uk/projects/fastqc
#
#	Documentation:	https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/
#
rule fastqc:
	input:
		unpack(fastqc_input)
	output:
		dir = directory(path.join(config['out'], 'fastqc/{sample}')),
		link_r1 = path.join(config['out'], 'fastqc/input', '{sample}_R1.fastq.gz'),
		link_r2 = path.join(config['out'], 'fastqc/input', '{sample}_R2.fastq.gz')
	log: set_log('fastqc/{sample}/{sample}_fastqc.log')
	shell:
		"""
		ln -s {input.read1} {output.link_r1}
		ln -s {input.read2} {output.link_r2}
		fastqc --noextract {output.link_r1} {output.link_r2} -o {output.dir} &> {log}
		"""

def multiqc_config(analysis_directory):
	"""
	Define custom multiqc options, create multiqc_config.yaml in analysis directory if it does not already exist

	:param analysis_directory:
	:return: full path to config file
	"""

	# multiqc configuration options
	multiqc_opts = [
		"fastqc_theoretical_gc: 'hg38_txome'"
	]
	multiqc_config_fp = os.path.join(analysis_directory, 'multiqc_config.yaml')

	# create multiqc_config yaml file if it doesn't already exist
	if not os.path.exists(multiqc_config_fp):
		cfg = 'fastqc_config:' + '\n  ' + '\n  '.join(multiqc_opts) + '\n'
		os.makedirs(analysis_directory, exist_ok=True)
		with open(multiqc_config_fp, 'w') as f:
			f.write(cfg)

	return multiqc_config_fp


# rule.fastqc
#
#	Description: 	Aggregate QC metrics into single report
#
# 	Reference: 		Ewels, P., Magnusson, M., Lundin, S., Käller, M., 2016. MultiQC: summarize analysis results for multiple tools
# 					and samples in a single report.
# 					Bioinformatics 32, 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
#
#	Documentation:	https://multiqc.info/docs
#
rule multiqc:
	input:
		multiqc_config(config['out']),
		expand(path.join(config['out'], 'fastqc/{sample}'), sample=SAMPLES),
		expand(
			path.join(config['out'], 'salmon/{sample}/quant.sf'),
			sample=SAMPLES
		) if opt_salmon else (),
		(
			expand(
				path.join(config['out'], 'star/{sample}/ReadsPerGene.out.tab'),
				sample=SAMPLES
			),
			expand(
				path.join(config['out'], 'star/{sample}/feature_counts.txt.summary'),
				sample=SAMPLES
			)

		) if opt_star else ()
	output:
		path.join(config['out'], "multiqc_report.html")
	params:
		'-m fastqc -m salmon -m star -m featureCounts --config ' + path.join(config['out'], 'multiqc_config.yaml')
	wrapper:
		"0.35.0/bio/multiqc"

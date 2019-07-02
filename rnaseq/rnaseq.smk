import os

from os import path
# todo: add validation step
from snakemake.io import protected, expand
from snakemake.io import directory
from snakemake import shell, rules
from snakemake.utils import R
from snakemake.utils import update_config, available_cpu_count

# r2py reports all stdout using warnings
# ignore warnings to maintain cleaner stdout
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)


# load analysis options
opt_star = config['star']
opt_salmon = config['salmon']

# set default configuration file
#configfile: 'config.yml'

SAMPLES = config['samples']
SAMPLE_IDS = list(SAMPLES.keys())

# location of salmon index directory
salmon_index_location =\
	os.path.join(
		os.path.dirname(config['transcripts_fa']),
		'salmon',
		os.path.splitext(os.path.basename(config['transcripts_fa']))[0]
	)

# todo: validate config input to ensure all files exist, generalize

# set default options
opts = dict()
opts['salmon-quant-l'] = 'ISR'
# override default options with config based options
if 'options' in config:
	opts.update(config['options'])


def set_log(name):
	return os.path.join(config['out'], name)

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

	#output:
	#	path.join(config['out'], 'dag.pdf')
	#shell: "snakemake -s rnaseq.smk --rulegraph | dot -Tpdf > {output}"

if opt_salmon:

	rule salmon_index:
		input:
			config['transcripts_fa']
		# index transcript fasta file if the directory does not already exist
		output:
			protected(
				directory(
					salmon_index_location
				)
			)
		log: 'salmon/salmon_index.log'
		shell:
			"""
			salmon index -t {input} -i {output} &> {log}
			"""

	# salmon libtype https://salmon.readthedocs.io/en/latest/salmon.html#what-s-this-libtype
	# pub: https://www.biorxiv.org/content/biorxiv/early/2016/08/30/021592.full.pdf
	rule salmon_quant:
		input:
			index = salmon_index_location,
			read1 = lambda w: config['samples'][w.sample]['read1'],
			read2 = lambda w: config['samples'][w.sample]['read2']
		output:
			quant = path.join(config['out'], 'salmon/{sample}/quant.sf')
		params:
			l = opts['salmon-quant-l']
		log: set_log('salmon/{sample}/stdout.log')

		run:
			output_dir = os.path.dirname(output[0])
			shell("salmon quant -i {input.index} -l {params.l} -1 {input.read1} -2 {input.read2} --validateMappings -o {output_dir} &> {log}")


	# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
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
# reference: https://github.com/alexdobin/STAR/blob/2.6.1d/doc/STARmanual.pdf
#
# optionally include STAR and downstream analysis components
#
if opt_star:

	## example config parameter values
	# genome_dir: /Precyte1/tmp/refs/genomes
	# genome_build: hg38
	# genome_id:  hg38wERCC92
	# genome_annotation_file: gencode.v25.primary_assembly.annotation.wERCC92.gtf
	# todo: check readlength, dynamically set sjdbOverhang
	# todo: allow mixed RL to be run in same analysis, currenly assumes the same for all samples
	# todo: create subdirectory based on GTF file, link or include GTF file with index
	rule star_index:
		input:
			fasta_files = path.join(config['genome_dir'], config['genome_build'], 'fa', config['genome_id']),
			annotation_gtf = config['genome_annotation_file']
		output:
			path = protected(
				directory(
					os.path.join(config['genome_dir'], config['genome_build'], 'indexes', config['genome_id'] + '_sjo' + str(config['star_sj_db_overhang']), 'star/')
				)
			)
		params:
			sj_db_overhang = str(config['star_sj_db_overhang'])
		# parameters:
		#	sjdbOverhang: 	'int>0: length of the donor/acceptor sequence on each side of the junctions'
		#		- value: set to ReadLength-1 (e.g. 75-1 for 75bp reads)
		#		- ref: STAR Manual (v2.6.1+)
		run:
			# check/create if output directory exists
			#if not os.path.exists(output.output_dir):
			#	os.makedirs(output.output_dir)
			command =\
			"star --runMode genomeGenerate --runThreadN 1 --genomeFastaFiles {input.fasta_files}/*.fa --genomeDir {output.path} --outFileNamePrefix {output.path}"
			# add optional parameters
			if input.annotation_gtf is not None and config['star_sj_db_overhang'] is not None:
				command += ' --sjdbGTFfile ' + input.annotation_gtf
				command += ' --sjdbOverhang ' + params.sj_db_overhang

			shell(command)

	rule star_align:
		input:
			read1 = lambda w: config['samples'][w.sample]['read1'],
			read2 = lambda w: config['samples'][w.sample]['read2'],
			#star_genome_dir = os.path.join(config['genome_dir'], config['genome_build'], 'indexes', config['genome_id'], 'star/')
			star_genome_dir = rules.star_index.output.path
		output:
			bam = path.join(config['out'], 'star/{sample}/Aligned.out.bam'),
			sorted_bam = path.join(config['out'], 'star/{sample}/Aligned.sortedByCoord.out.bam')
		params:
			quant_mode = 'TranscriptomeSAM GeneCounts',
			out_sam_type = 'BAM Unsorted SortedByCoordinate'
		threads: available_cpu_count()-2
		# parameters:
		#	quantMode: 'TranscriptomeSAM GeneCounts'
		#		- output: 	1) alignments translated into transcript coordinates
		#					2) number of reads/gene
		#		- ref: STAR Manual (v2.6.1+) - "Counting number of reads per gene."
		#	outSAMtype: 'BAM Unsorted SortedByCoordinate'
		#		- output: 	1) unsorted bam file
		#					2) sorted bam file
		#		- notes:	"Unsorted" output can be directly input into featureCounts
		shell:
			# option: --genomeLoad LoadAndKeep : osx incompatible
			"""
			star --runThreadN {threads} \
				 --genomeDir {input.star_genome_dir} \
				 --readFilesCommand gzcat \
				 --readFilesIn {input.read1} {input.read2} \
				 --outFileNamePrefix  $(dirname {output.sorted_bam})/ \
				 --quantMode {params.quant_mode} \
				 --outSAMtype {params.out_sam_type}
			"""


	rule feature_counts:
		input:
			#bam = path.join(config['out'], 'star/{sample}/Aligned.sortedByCoord.out.bam'),
			bam = rules.star_align.output.bam,
			annotation_gtf = config['genome_annotation_file']
		output:
			path.join(config['out'], 'star/{sample}/feature_counts.txt')
		params:
			strandedness = 2
		threads: available_cpu_count()-2
		# parameters:
		#	-p include if read are paired-end
		#	-B only count read-pairs that have both ends aligned
		#	-s 0|1|2 : unstranded|stranded|reversely stranded
		shell:
			"""
			featureCounts -T {threads} -p -B -s {params.strandedness} -a {input.annotation_gtf} -o {output} {input.bam}
			"""

# QC
rule fastqc:
	input:
		read1 = lambda w: config['samples'][w.sample]['read1'],
		read2 = lambda w: config['samples'][w.sample]['read2']
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

# see also https://multiqc.info/docs/#bulk-sample-renaming
rule multiqc:
	input:
		expand(path.join(config['out'], 'fastqc/{sample}'), sample=SAMPLES),
		expand(
			path.join(config['out'], 'salmon/{sample}/quant.sf'),
			sample=SAMPLES
		) if opt_salmon else ()
	output:
		path.join(config['out'], "multiqc_report.html")
	params:
		'-m fastqc -m salmon --config ' + path.join(config['out'], 'multiqc_config.yaml')
	wrapper:
		"0.35.0/bio/multiqc"

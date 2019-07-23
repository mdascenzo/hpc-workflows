

###### Cloning the repository:
```
git clone https://github.com/
```

###### Usage on Aloha:

The environment required to run workflows in this repository will be maintained on Aloha. No additional installation is needed. The example below 
shows the commands to activate the rnaseq environment on Aloha, and the basic commands to setup and execute the rnaseq workflow. A brief description 
of the configuration file and contents can be found later in this document. If using conda on Aloha for the first time, conda must be initialized.

One time initialization of conda. This will make an update to your bash_profile.
```
export PATH=/usr/local/env/conda/bin:$PATH
conda init
```

Activate rnaseq environment and execute workflow. The environment contains all the dependencies required by the workflow. 
```
# activate the analysis environment
conda activate aloha.rnaseq

# create an analysis directory
mkdir example_analysis

# download the rnaseq workflow within the analysis directory and link the workflow file (rnaseq.smk) and R directory 
cd example_analysis
git clone https://github.com/
ln -s workflows/src/rnaseq/* .

# create an analysis configuration file using compressed paired-end fastq reads as input
create_config.py -f path/to/sequence/*fq.gz -o analysis

# examine the newly created config.yml file to customize the analysis (see below). 

# execute workflow using 8 cores
snakemake -s rnaseq.smk --cores 8
```

###### Conda Installation on OSX or Linux

The commands below can be used to create a new conda environment on either OSX or Linux.
```
# set install directory (e.g. /usr/local/env/conda) and set as current directory
INSTALL_DIR=/usr/local/env/conda
cd $INSTALL_DIR

# for linux
curl https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o Miniconda3-3.7.0-Linux-x86_64.sh

# for osx
curl https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o Miniconda3-latest-MacOSX-x86_64.sh

# execute installer including desire install location (-p PREFIX):
sh Miniconda3-3.7.0-Linux-x86_64.sh -p $INSTALL_PATH
 
# update path to include the newly created conda install location (e.g. /usr/local/env/conda) 
export PATH=$INSTALL_DIR/bin:$PATH

# install anaconda client
conda install anaconda-client
```

Once Conda is installed, all required dependencies can be installed using the following commands.
```
git clone https://github.com/
cd workflows 
conda env create -f env/rnaseq.yml
```

###### To use an alternate configuration file:
```
snakemake --configfile <filename.yml> -s rnaseq.smk --cores 8
```

###### Updating an existing installation:
```
conda env update --name rnaseq  -f env/rnaseq.yml
```

###### Example config.yml file:

The configuration file contains default settings and should be customized. The create_config.py setup script includes all paired-end fastq.gz files
under the "samples" sub-section of the YAML document. 
```
analysis_name: rnaseq_analysis

resources_dir: /Volumes/Precyte1/stage/resources
build: hg38
genome_uid: hg38wERCC92

star: yes
star_sj_db_overhang: 74
annotation_gtf: gencode.v25.primary_assembly.annotation.wERCC92.gtf

salmon: yes
tx_uid: ensembl_rel83
tx2gene_fp: /Precyte1/stage/resources/transcriptomes/hg38/ensembl_rel86/annotation/tx2gene/tx2gene.EnsDb.Hsapiens.v86.csv

trim: yes

out: /example/path/DATE_example/analysis

options:
  trimmomatic-adapters-fa: TruSeq3-PE-2.fa

samples:
  BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R1_001:
    read1: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R1_001.fastq.gz
    read2: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R2_001.fastq.gz
  BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R1_001:
    read1: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R1_001.fastq.gz
    read2: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R2_001.fastq.gz
```

Note that configuration settings define paths to resources currently located in /Precyte1/stage/resources.

###### Analysis Workflow Options
Example workflow options (v0.0.1). Trim step not currently shown.

![workflow-full](doc/rnaseq/img/dag_options_3.png)

#### Notes
=
##### STAR

Genome Indexing:
- Configured to utilize a splice junction database.
- sjdbOverhang option should be set to ReadLength-1 (e.g. 74 for 75bp PE reads)

options:
- sjdbGTFfile
- sjdbOverhang

Example splice junction db file:
```
/Precyte1/stage/refs/genomes/hg38/annotation/hg38wERCC92/gencode.v25.primary_assembly.annotation.wERCC92.gtf
```

###### Todo:
- Possibly create installer to move R-code to library path
- Update genome and transcriptome annotation files

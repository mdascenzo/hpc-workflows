

###### Cloning the repository:
```
git clone https://github.com/ 
```

###### Installation

If installing for the first time:
```
cd workflows
conda env create -f env/rnaseq.yml
```

###### Example usage:
```
# activate rnaseq environment
conda activate rnaseq

# create analysis config file
create_config.py -f seq/*fq.gz -o analysis

# execute workflow
snakemake -s rnaseq.smk --cores 8
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

out: /example/path/DATE_example/analysis

options:
  trimmomatic-adapters-fa: TruSeq3-PE-2.fa

samples:
  BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R1_001:
    read1: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R1_001.fastq.gz
    read2: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044190_ATTACTCG-TATAGCCT_L006_R2_001.fastq.gz
    trim: no
  BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R1_001:
    read1: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R1_001.fastq.gz
    read2: /example/path/DATE_example/seq/BHJ3VCBBXX_WB-1044191_ATTACTCG-ATAGAGGC_L006_R2_001.fastq.gz
    trim: no
```

###### Analysis Workflow Options
![workflow-full](doc/rnaseq/img/dag_options_3.png)

###### Todo:
- Possibly create installer to move R-code to library path
- Update genome and transcriptome annotation files

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
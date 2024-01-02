
This repository contains utilities for implementing HPC (High-Performance Computing) based workflows on [AWS ParallelCluster](https://aws.amazon.com/hpc/parallelcluster/) in conjunction with the [Slurm](https://slurm.schedmd.com/documentation.html) workload manager. A [Snakemake](https://snakemake.readthedocs.io/en/stable/)-based RNA-Seq workflow is included, along with example [Packer](https://www.packer.io) build scripts for automating the creation of a custom AMI. This AMI contains the required RNA-Seq analysis software, bundled with the necessary dependencies to run within the AWS ParallelCluster environment.

## Workflows
### RNA-Seq

![workflow-full](doc/rnaseq/img/dag_options_3.png)

## Analysis Configuration

## Example config.yml:
```yml
analysis_name: rnaseq_analysis
star: yes
salmon: yes
trim: yes
resources_dir: /nfs/resources
build: hg38
genome_uid: hg38wERCC92
tx_uid: ensembl_rel83
annotation_gtf: gencode.v25.primary_assembly.annotation.wERCC92.gtf
tx2gene_fp: /nfs/resources/transcriptomes/hg38/ensembl_rel86/annotation/tx2gene/tx2gene.EnsDb.Hsapiens.v86.csv
star_sj_db_overhang: 149
out: /nfs/workspace/out/run001/analysis
options:
  trimmomatic-adapters-fa: TruSeq3-PE-2.fa
samples:
  S001:
    read1: /nfs/workspace/run001/S001_R1.fastq.gz
    read2: /nfs/workspace/run001/S001_R2.fastq.gz
  S002:
    read1: /nfs/workspace/run001/S002_R1.fastq.gz
    read2: /nfs/workspace/run001/S002_R2.fastq.gz
```
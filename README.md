# Introduction

RNA-Seq workflow optimized for use with [AWS Parallel Cluster](https://aws.amazon.com/hpc/parallelcluster/) and [Slurm workload manager](https://slurm.schedmd.com/documentation.html).

This branch adds functionality that integrates RNA-Seq workflow with AWS Parallel Cluster and Slurm workload manager. Use of Slurm allows greater control over compute resources and creates coherent log files. The AWS Parallel Cluster configuration included in this branch utilizes spot instances offering a significant cost reduction compared to on-demand instanses. Additionally, compute nodes are autoscaled to meet analysis workload and are terminated when not in use for 10 minutes. This branch includes [Packer](https://www.packer.io) build scripts for automating the build of a custom AMI containing the required RNA-Seq analysis dependencies/tools alongside AWS Parallel Cluster and Slurm resources. Docker is no longer implemented in this branch but may be re-implemented in future versions.

This is a development branch. Optimizations to workflow resources (memory, cpus, parallel tasks per node, etc.) will likely be made as this workflow is applied to new data.

# Notes for use on AWS:

## Prerequisites:
- Python 3.0 
- AWS ParallelCluster installed on a local machine (see: [Installing AWS ParallelCluster](https://docs.aws.amazon.com/parallelcluster/latest/ug/install.html))
- 

## Quickstart
After installing prerequisites, copy the pcluster configuration to your home directory:
```
git clone https://github.com/
git checkout dev-hpc
cp workflows/aws/parallelcluster/config ~/.parallelcluster
```

#### Create a new cluster
```
pcluster create rnaseq
```
This command takes about 5-10 minutes to complete as it provisions a Master node, Compute nodes, and EBS backed NFS resources for sharing data between nodes. By default Compute nodes are not launched, but will launch once the workflow is started. This allows data to be copied from S3 to the cluster, avoiding idle Compute nodes. 

#### Connect to the head node
```
pcluster ssh rnaseq
```
NFS mounts are auto created:
- /workspace
- /research

#### Copy data from S3 to cluster
```
aws s3 cp --recursive s3:// /workspace
```
Resource files (STAR index files, annotation, etc) are auto mounted in /research and do not need to be copied to the cluster. 

### Install workflow setup script

Workflow setup scripts are not currently installed by default on the cluster, install on the head node:
```
git clone https://github.com/ 
git checkout dev-hpc
cd workflows
sudo pip install .
```

#### Run RNA-Seq Workflow
```
# create and cd to analysis/working directory
analysis_dir='/workspace/30-410354445'
cd $analysis_dir
mkdir seq
mv *.fastq.gz ./seq

Workflow setup scripts are not currently installed by default on the cluster, install now:
git clone https://github.com/ 
git checkout dev-hpc
cd workflows
sudo pip install .

# link contents of the rnaseq working directory to the current analysis/working dir
# contents include the rnaseq.smk workflow and required R scripts
cd $analysis_dir
ln -s workflows/src/rnaseq/* .

# create analysis configuration file
# the default configuration needs to be updated, a working example is below.
create_config.py -f seq/*.gz

# run the analysis using Slurm profile
nohup snakemake -s rnaseq.smk --profile slurm &> out.log &

# to gracefully stop workflow:
killall -TERM snakemake
```

#### Example config.yml header:
```yml
analysis_name: rnaseq_analysis
star: yes
salmon: yes
trim: yes
resources_dir: /research/resources
build: hg38
genome_uid: hg38wERCC92
tx_uid: ensembl_rel83
annotation_gtf: gencode.v25.primary_assembly.annotation.wERCC92.gtf
tx2gene_fp: /research/resources/transcriptomes/hg38/ensembl_rel86/annotation/tx2gene/tx2gene.EnsDb.Hsapiens.v86.csv
star_sj_db_overhang: 149
out: /workspace/30-410354445/analysis
options:
  trimmomatic-adapters-fa: TruSeq3-PE-2.fa
samples:
  01-08-14-20_R1_001:
    read1: /workspace/30-410354445/seq/01-08-14-20_R1_001.fastq.gz
    read2: /workspace/30-410354445/seq/01-08-14-20_R2_001.fastq.gz
```
# Introduction

RNA-Seq workflow optimized for use with [AWS Parallel Cluster](https://aws.amazon.com/hpc/parallelcluster/) and [Slurm workload manager](https://slurm.schedmd.com/documentation.html).

This branch adds functionality that integrates RNA-Seq workflow with AWS Parallel Cluster and Slurm workload manager. Use of Slurm allows greater control over compute resources and creates coherent log files. The AWS Parallel Cluster configuration included in this branch utilizes spot instances offering a significant cost reduction compared to on-demand instanses. Additionally, compute nodes are autoscaled to meet analysis workload and are terminated when not in use for 10 minutes. This branch includes [Packer](https://www.packer.io) build scripts for automating the build of a custom AMI containing the required RNA-Seq analysis dependencies/tools alongside AWS Parallel Cluster and Slurm resources. 

Docker is implemented to provide an environement for launching and managing the HPC cluster. 

This is a development branch. Optimizations to workflow resources (memory, cpus, parallel tasks per node, etc.) will likely be made as this workflow is applied to new data.

### Prerequisites:

The following files must be in place and configured prior to proceeding to subsequent steps.

Setting | Location | Path
--- | --- | ---
AWS config | local |  ~/.aws/config
AWS credentials | local |  ~/.aws/credentials
ssh config | local | ~/.ssh/config
ssh credentials | local |  ~/.ssh/

These files are are mounted within 

##Cluster Setup

### Local Docker Environment
Docker is used to maintain the required environment for manage thecluster. Create a Docker container using the following image:
```
docker pull X
docker run --entrypoint /bin/bash
```

TODO: copy to docker container
git clone https://github.com/ 
git checkout dev-hpc                                          
cp workflows/aws/parallelcluster/config ~/.parallelcluster 

#### Create cluster
From within the docker container:
```
pcluster create
```
This command takes about 5-10 minutes to complete as it provisions a Master node, Compute nodes, and EBS backed NFS resources for sharing data between nodes. By default Compute nodes are not launched, but will launch once the workflow is started. This allows data to be copied from S3 to the cluster, avoiding idle Compute nodes. 

#### Configure SSH timeout
To avoid timeouts when connection to the head node via SSH, update your SSH config (located at ~/.ssh/config) to contain the following setting:

`ServerAliveInterval 120`

#### ADD SSH credentials to head node
Copy credentials to head node for use with aws command line tools (awscli).
```
HEAD_NODE_IP_ADDRESS=`aws ec2 describe-instances --query "Reservations[*].Instances[*].PublicIpAddress" --filter Name=tag:Name,Values=Master --output text`
scp -r -i ~/.ssh/ ~/.aws ubuntu@${HEAD_NODE_IP_ADDRESS}:.aws
```

### Copy data
Copy sequence data and analysis resource files using awscli.
```
# List the contents of the s3 shared folder:
aws s3 ls s3://

# Use the desired folder name to set the value of ANALYSIS_ID, for example:
ANALYSIS_ID=30-410354445

# Sync the folder set by ANALYSIS_ID in the previous command to EC2 workspace
aws s3 sync s3:// /workspace/${ANALYSIS_ID}
```

Copy resource files (STAR index files, annotation, etc) to /research:
```
aws s3 sync s3:// /research/resources
```

#### Run RNA-Seq Workflow
```
# create and cd to analysis/working directory
ANALYSIS_DIR='/workspace/${ANALYSIS_ID}'
cd $ANALYSIS_DIR
mkdir seq
mv *.fastq.gz ./seq

# workflow setup scripts are not currently installed by default on the cluster, install now:
git clone https://github.com/ 
git checkout dev-hpc
cd workflows
sudo pip install .

# link contents of the rnaseq working directory to the current analysis/working dir
# contents include the rnaseq.smk workflow and required R scripts
cd $ANALYSIS_DIR
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
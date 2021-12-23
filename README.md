# Introduction

RNA-Seq workflow optimized for use with [AWS Parallel Cluster](https://aws.amazon.com/hpc/parallelcluster/) and [Slurm workload manager](https://slurm.schedmd.com/documentation.html).

This branch adds functionality that integrates RNA-Seq workflow with AWS Parallel Cluster and Slurm workload manager. Use of Slurm allows greater control over compute resources and creates coherent log files. The AWS Parallel Cluster configuration included in this branch utilizes spot instances offering a significant cost reduction compared to on-demand instanses. Additionally, compute nodes are autoscaled to meet analysis workload and are terminated when not in use for 10 minutes. This branch includes [Packer](https://www.packer.io) build scripts for automating the build of a custom AMI containing the required RNA-Seq analysis dependencies/tools alongside AWS Parallel Cluster and Slurm resources. 

Docker is implemented to provide an environement for launching and managing the HPC cluster. 

This is a development branch. Optimizations to workflow resources (memory, cpus, parallel tasks per node, etc.) will likely be made as this workflow is applied to new data.

## Initial Setup

In the following steps a Docker environement will be created to facilate creation and management of an HPC cluster used to process analysis workflows. AWS and ssh credentials must be setup to allow communication between a local machine (e.g. Aloha or iMac) and AWS.

### AWS User Permissions

Users should first be assigned to the **IAM User** group `PrecyteClusterUsers` within the AWS Control Panel.

### Repository
Clone the workflow repository to a local machine.
```
git clone https://github.com/

# - or - 

git clone git@github.com:

# checkout the hpc branch:
git checkout dev-hpc
```

### Configuration files
The following configuration files are required before proceeding. The files should be created on a **local machine** within your home directory, later the parent directoreis (.aws, .parallelcluster, and .ssh) will be mounted inside a Docker container launched in subsequent steps.

Setting | Location | Path
--- | --- | ---
AWS config | local |  ~/.aws/config
AWS credentials | local |  ~/.aws/credentials
AWS parallel cluster config | local | ~/.parallelcluster/config 
ssh config | local | ~/.ssh/config
ssh credentials | local |  ~/.ssh/

### AWS Setup
On local machine (e.g. iMac, Aloha):

#### AWS config setup (~/.aws/config):

```
[default]
region = us-west-2
output = json
```

#### AWS credentials setup (~/.aws/credentials):
```
[default]
aws_access_key_id = AKIARGXPHLFIQCXXXXXX
aws_secret_access_key = u2p16+VicWwhqFRyaTbeZO5SVbOgBOtDscXXXXXX
```

#### AWS Parallel Cluster setup (~/.parallelcluster/config):
```
mkdir ~/.parallelcluster

# from within workflows directory
cp aws/parallelcluster/config ~/.parallelcluster 
```

### SSH Setup

On local machine (e.g. iMac, Aloha):

#### .ssh directory setup

```
# create .ssh directory if needed and set permissions
mkdir ~/.ssh
chmod 700 ~/.ssh

# copy to local .ssh directory and set permissions
cp /Volumes/Precyte1/research/aws/ssh/ ~/.ssh
chmod 400 ~/.ssh/
```

#### SSH config setup (~/.ssh/config):
To avoid timeouts when connecting to the head node via SSH, update your SSH config to contain the following setting:

`ServerAliveInterval 120`

The following command will add this setting and create the config file if it is not already present.
```
echo "ServerAliveInterval 120" >> ~/.ssh/config 
```

### Local Docker Setup

Docker is used to maintain the required environment for managing thecluster Before proceeding, check that the following are installed on your local machine. For Docker on OSX, Docker Compose V2 must be enabled within the Docker control panel.

#### Requirements:
- Docker
- Docker Compose (version v2.2.1 or greater) 

Verify Docker Compose version with: `docker-compose --version`

## Launch Cluster

To do: add a description of the multiple layers utilized to access AWS.

The layers are:

1. Launch environment: A Docker container containing software required to manage cluster via AWS Parallel Cluster using pcluster. Provies utility to start, stop, delete cluster and acts as a gateway to connect to the cluster head node.
2. HPC cluster 

### Create launch environment:

```
# from within the workflow repository
cd workflows

# build and start docker container
docker-compose -f .docker/docker-compose.yaml up

# in a new terminal, connect to the container
docker exec -it pcluster_admin /bin/bash
```

### Create cluster
From within the launch environment (docker container):
```
pcluster create
```
This command takes about 5-10 minutes to complete as it provisions a Head node, Compute nodes, and EBS backed NFS resources for sharing data between nodes. By default, Compute nodes are not launched, but will launch once the workflow is started allowing data to be copied from S3 to the cluster while avoiding idle Compute nodes. 

#### Copy AWS credentials and workflows to head node:
Once the cluster has started, copy the following to the head node:
 1. AWS credentials, for use with AWS command line tools (awscli)
 2. /code/workflows

The path `/code/workflows` should have been auto-mounted for quick access within the launch environment. 

```
HEAD_NODE_IP_ADDRESS=`aws ec2 describe-instances --query "Reservations[*].Instances[*].PublicIpAddress" --filter Name=tag:Name,Values=Master --output text`

# copy .aws directory to head node
scp -r -i ~/.ssh/ ~/.aws ubuntu@${HEAD_NODE_IP_ADDRESS}:.aws

# copy workflow repository to head node, exclude hidden files
rsync -a --progress  -e "ssh -i ~/.ssh/ --exclude='.*'  /code/workflows ubuntu@${HEAD_NODE_IP_ADDRESS}:/workspace
```
## Setup and Run Analysis
### SSH to head node
Login to the head node from the launch environment using:
```
pcluster ssh
```

### Install workflows
Analysis workflow setup scripts were copied to the cluster in a previous step. These are not currently installed by default on the cluster, install now:
```
cd /workspace/workflows
sudo pip install .
pip install --user .
export PATH=$PATH:/home/ubuntu/.local/bin
```

### Copy data
#### Copy sequence data and analysis resource files using awscli.
```
# List the contents of the s3 shared folder:
aws s3 ls s3://

# Use the desired folder name to set the value of ANALYSIS_ID, for example:
ANALYSIS_ID=30-410354445

# Sync the folder set by ANALYSIS_ID in the previous command to EC2 workspace
aws s3 sync s3:// /workspace/${ANALYSIS_ID}
```

#### Copy resource files (STAR index files, annotation, etc) to /research:
```
aws s3 sync s3:// /research/resources
```

### Run RNA-Seq Workflow
### Run RNA-Seq Workflow
```
# create and cd to analysis/working directory
ANALYSIS_DIR=/workspace/${ANALYSIS_ID}
cd $ANALYSIS_DIR
mv 00_fastq seq

# link contents of the rnaseq working directory to the current analysis/working dir
# contents include the rnaseq.smk workflow and required R scripts
cd $ANALYSIS_DIR
ln -s /workspace/workflows/src/rnaseq/* .

# create analysis configuration file
# the default configuration needs to be updated, a working example is below.
create_config.py -f seq/*.gz -o analysis

# run the analysis using Slurm profile
nohup snakemake --keep-going -s rnaseq.smk --profile slurm &> out.log &
```

### Monitor analysis progress
```
grep done$ /workspace/30-601341635/out.log | tail -n 1
```


### Copy completed analysis data to S3
```
# sync to
aws s3 sync ${ANALYSIS_DIR} s3:// \
  --exclude 'bam/*' \
  --exclude 'seq/*' \
  --exclude '*trim_reads/*' \
  --exclude "*.snakemake/*" --exclude "*/.git/*" \
  --exclude "*.bam*" --exclude "*.fq" --exclude "*.fastq" --exclude "*.fastq.gz" --exclude "*.fq.gz" 
  ```

---

## Useful Commands

### snakemake:

#### killall 
Gracefully stop workflow. All queued jobs will run to completion.
```
killall -TERM snakemake
```
---
---

### slurm commands:

#### sinfo
Display available compute resources including CPU load:
```
sinfo -o "%20P %.5a %.10l %.6D %.6t %40N  %5c %20C %20O "
```

#### squeue
Display queue information with added width showing full job name:
```
squeue -o "%.18i %.9P %.50j %.8u %.2t %.10M %.6D %R"
```

See: [multiple-queue-mode-slurm-user-guide](https://docs.aws.amazon.com/parallelcluster/latest/ug/multiple-queue-mode-slurm-user-guide.html)

#### scontrol
Show additional job details:
```
scontrol show jobid -dd <job-id>
```
--- 

## Other Commands

Check analysis progress:
```
 grep done$ /workspace/30-601341635/out.log | tail -n 1
```

## Admin Commands

Log into compute node from head node.
```
ssh compute1-dy-c5a8xlarge-1
```

Cancel all jobs
```
squeue --me -h -o "%i" | xargs scancel
```

## Change cluster status
For example, change from spot to ondemand. Update config file in ~/.parallelcluster/config, then:
```
pcluster update
```



## Example config.yml header:
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


## Troubleshooting

#### Jobs submitted to cluster but compute nodes do not start.
Spot instances may not be available, confirm using [cloudfront_dashbord:slurm_resume](https://us-west-2.console.aws.amazon.com/cloudwatch/home?region=us-west-2#dashboards:name=parallelcluster-

Verify the following error:
ERROR - Encountered exception when launching instances for nodes (x2) ['compute1-dy-c5a8xlarge-1', 'compute1-dy-c5a8xlarge-2']: 
An error occurred (InsufficientInstanceCapacity) when calling the RunInstances operation (reached max retries: 1): There is no Spot capacity available that matches your request.


## Notes
```
sudo systemctl stop unattended-upgrades
sudo apt-get purge unattended-upgrades
```

#### Transfer speed
```
# from head node
dd if=/dev/zero of=/workspace/tmp.img bs=1G count=1 oflag=direct
1+0 records in
1+0 records out
1073741824 bytes (1.1 GB, 1.0 GiB) copied, 45.0649 s, 23.8 MB/s
```

From compute node
```
1073741824 bytes (1.1 GB, 1.0 GiB) copied, 37.9573 s, 28.3 MB/s
```

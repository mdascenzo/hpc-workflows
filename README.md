
# Notes for use on AWS:

## Instance:

Tested on r5a instance type based on AMD EPYC 7000. [[link]](https://aws.amazon.com/ec2/instance-types/r5/) 

* Minimum required resources: 16 vCPU / 128 Gb RAM | **r5a.4xlarge**
* For larger analyses, consider: 48 vCPU / 384 Gb RAM | **r5a.12xlarge**

A template has been created and configured to run instances using a spot request. This has pros/cons, on the upside pricing is favorable (e.g r5ad.12xlarge @ $0.85 spot vs. $3.14 on-demand), however availability can sometimes be limited, which can result in termination of the instance mid-analysis. While termination is less than ideal, the analysis is backed by a persistent EBS volume, allowing the analysis to be resumed after re-launching a new instance.

* Template [[link]](https://us-west-2.console.aws.amazon.com/ec2/v2/home?region=us-west-2#LaunchTemplateDetails:launchTemplateId=lt-09c3433e9b361f599)

## Storage:

EBS attached storage is used for analysis resource (e.g. STAR index and annotation), FASTQ input files, and analysis output. 
* Size: As a general guideline, 1500 GB should be adequate for 25 samples.
* Recommended storage type: Throughput Optimized (St1)


#### Attach volume:
Newly created volume must be manually attached and mounted:
```
# check storage name
lsblk

# format
sudo mkfs -t xfs /dev/nvme1n1

# create mount point
sudo mkdir /nfs

# using nvme1n1 based on lsblk, however this is dynamic
sudo mount -t xfs /dev/nvme1n1 /nfs
```

## AWS Credentials

#### Add credentials
```
aws configure

AWS Access Key ID [None]: AKIAIOSFODNN7EXAMPLE
AWS Secret Access Key [None]: wJalrXUtnFEMI/K7MDENG/bPxRfiCYEXAMPLEKEY 
Default region name [None]: us-west-2 
Default output format [None]: json
```

### Monitor cost of spot instance

```
#start:
aws ec2 create-spot-datafeed-subscription --bucket --prefix spot/

#stop:
aws ec2 delete-spot-datafeed-subscription
```

## Docker

Docker is already setup and configured in the launched AMI. It should contain the most recent rnaseq Docker image. If a newer version is available, it can be pulled from Amazon Elastic Container Registry (ECR) using the following steps, otherwise, skep to step 4 to launch the existing docker image installed in the launched AMI. 

#### Pull Docker Image:
```
docker pull.amazonaws.com/rnaseq:0.1.7
```

#### Create Docker Container:

This will create as a daemon which will allow the container to continue running even after existing the shell. 

```
docker run -dt -v /nfs:/nfs -v /home/admin/.aws:/root/.aws --name rnaseq.amazonaws.com/rnaseq:0.1.7
```

To check that the container is running. 
```
docker ps
```

To connect to the container, this may take a bit (~30s) to connect:
```
docker exec -it rnaseq /bin/bash 
```

## S3 / Data

Sequence files must be manually copied from S3. For example:

```
aws s3 cp --recursive s3:// .
```
!! Note: Resources files must also be copied over from S3 (STAR index files, annotation, etc).

## RNA-Seq Workflow

All commands executed to run the RNA-Seq workflow should be run from within the launched rnaseq Docker container. 

From within docker container:
```
# create and change analysis/working directory
mkdir analysis_dir
cd analysis_dir

# clone workflow from GitHub:
git clone https://github.com/ 
git checkout remotes/origin/dev-docker

# link contents of the rnaseq working directory to the current analysis/working dir
# contents include the rnaseq.smk workflow and required R scripts
ln -s workflows/src/rnaseq/* .

# create analysis configuration file
# the default configuration needs to be updated, a working example is below.
create_config.py -f seq/*.gz

# run the analysis 
nohup snakemake -s rnaseq.smk &> out.log &

# to gracefully stop workflow:
killall -TERM snakemake
```

#### Example config.yml header:
```yml
analysis_name: rnaseq_analysis
star: yes
salmon: yes
trim: yes
resources_dir: /nfs/research/resources
build: hg38
genome_uid: hg38wERCC92
tx_uid: ensembl_rel83
annotation_gtf: gencode.v25.primary_assembly.annotation.wERCC92.gtf
tx2gene_fp: /nfs/research/resources/transcriptomes/hg38/ensembl_rel86/annotation/tx2gene/tx2gene.EnsDb.Hsapiens.v86.csv
star_sj_db_overhang: 149
out: /nfs/30-410354445/analysis
options:
  trimmomatic-adapters-fa: TruSeq3-PE-2.fa
samples:
  01-08-14-20_R1_001:
    read1: /nfs/30-410354445/seq/01-08-14-20_R1_001.fastq.gz
    read2: /nfs/30-410354445/seq/01-08-14-20_R2_001.fastq.gz
```
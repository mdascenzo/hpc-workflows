#!/bin/bash

## set environment variables

# set path for install
export DEBIAN_FRONTEND=noninteractive
#export PATH="/usr/local/env/conda/bin:${PATH}"
# set path for user
#sudo echo "PATH=\"/usr/local/env/conda/bin:${PATH}\"" > /etc/environment
##
#exit 0
apt-get update
# remove upgrades that autorun daily
apt-get remove -y unattended-upgrades
# turn off motd
chmod -x /etc/update-motd.d/*

# build tools, python, and R, sysstat
apt-get install -y \
    build-essential cmake git swig libz-dev libbz2-dev vim \
    tree htop wget zip rsync pigz rename sysstat \
    graphviz libgraphviz-dev pkg-config 

# java
apt-get install -y default-jdk

# python (conda)
# see: https://docs.conda.io/en/latest/miniconda.html
mkdir /usr/local/sw
mkdir /tmp/build
cd /tmp/build
wget https://repo.anaconda.com/miniconda/Miniconda3-py311_23.11.0-2-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p /usr/local/env/conda

# export paths for this build
export PATH="/usr/local/env/conda/bin:${PATH}"

# install mamba
conda init
conda config --add channels anaconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install -y -c conda-forge mamba

# general packages
#conda install -y -c anaconda pandas #graphviz

# python packages
mamba install -y pandas=2.1.4 cython=3.0.7

# R and required packages
mamba install -y -c conda-forge r-base=4.3.2 r-tidyverse=2.0.0 r-hmisc=5.1_1 r-kernsmooth=2.23_22 r-rlist=0.4.6.2 r-locfit=1.5_9.8
R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org', quiet=FALSE, verbose=TRUE)"
R -e "BiocManager::install(version = '3.18')"
R -e "BiocManager::install(c('DESeq2', 'tximport', 'dupRadar'))"
mamba install -y -c conda-forge rpy2=3.5.11

# STAR
cd /tmp/build
wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz
tar -xzf 2.7.3a.tar.gz
cd STAR-2.7.3a
cd source
make STAR
mv STAR /usr/local/bin

# snakemake
cd /tmp/build
git clone https://github.com/snakemake/snakemake.git
cd snakemake
git checkout v5.10.0
python setup.py install
ln -s /usr/local/env/conda/bin/snakemake /usr/local/bin/sm

# salmon
cd /tmp/build
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz
tar -xvf salmon-1.1.0_linux_x86_64.tar.gz
mv salmon-latest_linux_x86_64 /usr/local/sw/salmon-1.1.0
ln -s /usr/local/sw/salmon-1.1.0/bin/salmon /usr/local/bin

# trimmomatic
cd /tmp/build
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip Trimmomatic-0.39.zip
mv Trimmomatic-0.39 /usr/local/sw
ln -s /usr/local/sw/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin

# subread package (only featureCounts is linked)
cd /tmp/build
wget --tries=70 https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz
tar -xvf subread-2.0.0-Linux-x86_64.tar.gz
mv subread-2.0.0-Linux-x86_64 /usr/local/sw/subread-2.0.0
ln -s /usr/local/sw/subread-2.0.0/bin/featureCounts /usr/local/bin

# spades
cd /tmp/build
wget http://cab.spbu.ru/files/release3.14.0/SPAdes-3.14.0-Linux.tar.gz
tar -xvf SPAdes-3.14.0-Linux.tar.gz
mv SPAdes-3.14.0-Linux /usr/local/sw/SPAdes-3.14.0
ln -s /usr/local/sw/SPAdes-3.14.0/bin/* /usr/local/bin

# FastQC
cd /tmp/build
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
chmod 775  FastQC/fastqc
mv FastQC /usr/local/sw
ln -s /usr/local/sw/FastQC/fastqc /usr/local/bin

# multiqc
cd /tmp/build
git clone https://github.com/ewels/MultiQC.git
cd MultiQC
git checkout v1.8
python setup.py install

# picard
cd /tmp/build
wget https://github.com/broadinstitute/picard/releases/download/2.24.1/picard.jar
mv picard.jar /usr/local/sw
ln -s /usr/local/sw/picard.jar /usr/local/bin

# samtools
apt-get -y install libssl-dev libcurl4-gnutls-dev zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev
cd /tmp/build 
wget https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 
tar -vxjf htslib-1.11.tar.bz2 
cd htslib-1.11 
./configure --prefix /usr/local/sw/htslib-1.11 
make 
make install 
ln -s /usr/local/sw/htslib-1.11/bin/* /usr/local/bin 
cd /tmp/build 
wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 
tar -vxjf samtools-1.11.tar.bz2  
cd samtools-1.11 
./configure --prefix /usr/local/sw/samtools-1.11 
make 
make install 
ln -s /usr/local/sw/samtools-1.11/bin/* /usr/local/bin 
cd /tmp/build 
wget https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 
tar -vxjf bcftools-1.11.tar.bz2 
cd bcftools-1.11 
./configure --prefix /usr/local/sw/bcftools-1.11 
make 
make install 
ln -s /usr/local/sw/bcftools-1.11/bin/* /usr/local/bin
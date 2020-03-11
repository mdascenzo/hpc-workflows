FROM debian:latest

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

# install build tools, python, and R
RUN apt-get install -y \
  build-essential cmake git python-pip swig libz-dev libbz2-dev vim \
  tree htop wget zip rsync pigz rename

# install java
RUN apt-get install -y default-jdk

# install python (conda)
RUN mkdir /usr/local/sw && \
    mkdir /tmp/build && \
    cd /tmp/build && \
    wget https://repo.continuum.io/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh -O miniconda.sh && \
    sh miniconda.sh -b -p /usr/local/env/conda

# setup paths
ENV PATH="/usr/local/env/conda/bin:${PATH}"

# update bashrc
RUN echo "alias ll='ls -l'" >>  ~/.bashrc

# install anconda
RUN conda init && \
    conda install -y anaconda-client

# install python packages
RUN conda install -y pandas

# R 3.5.2 and required libraries
#RUN apt-get qinstall -y r-base=3.5.2-1 && \
#    pip install rpy2==3.2.6 && \
#    R -e "install.packages('tidyverse', repos = 'https://cloud.r-project.org', quiet=FALSE, verbose=TRUE)" # non-versioned

RUN conda install -y -c r r-base=3.6.1 r-tidyverse=1.2.1&& \
    conda install -y -c conda-forge r-rlist=0.4.6.1 && \
    R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org', quiet=FALSE, verbose=TRUE)" && \
    R -e "BiocManager::install(c('DESeq2', 'tximport', 'dupRadar'))"  && \
    conda install -y rpy2 && \
    conda install -y -c r r-kernsmooth

# STAR
RUN cd /tmp/build && \
    wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz && \
    tar -xzf 2.7.3a.tar.gz && \
    cd STAR-2.7.3a && \
    cd source && \
    make STAR && \
    mv STAR /usr/local/bin

# snakemake
RUN cd /tmp/build && \
    git clone https://github.com/snakemake/snakemake.git && \
    cd snakemake && git checkout v5.10.0 && \
    python setup.py install && \
    ln -s /usr/local/env/conda/bin/snakemake /usr/local/bin/sm

# salmon
RUN cd /tmp/build && \
    wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz && \
    tar -xvf salmon-1.1.0_linux_x86_64.tar.gz && \
    mv salmon-latest_linux_x86_64 /usr/local/sw/salmon-1.1.0 && \
    ln -s /usr/local/sw/salmon-1.1.0/bin/salmon /usr/local/bin

# trimmomatic
RUN cd /tmp/build && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv Trimmomatic-0.39 /usr/local/sw && \
    ln -s /usr/local/sw/Trimmomatic-0.39/trimmomatic-0.39.jar /usr/local/bin

# subread package (only featureCounts is linked)
RUN cd /tmp/build && \
    wget https://sourceforge.net/projects/subread/files/subread-2.0.0/subread-2.0.0-Linux-x86_64.tar.gz && \
    tar -xvf subread-2.0.0-Linux-x86_64.tar.gz && \
    mv subread-2.0.0-Linux-x86_64 /usr/local/sw/subread-2.0.0 && \
    ln -s /usr/local/sw/subread-2.0.0/bin/featureCounts /usr/local/bin

# spades
RUN cd /tmp/build && \
    wget http://cab.spbu.ru/files/release3.14.0/SPAdes-3.14.0-Linux.tar.gz && \
    tar -xvf SPAdes-3.14.0-Linux.tar.gz && \
    mv SPAdes-3.14.0-Linux /usr/local/sw/SPAdes-3.14.0 && \
    ln -s /usr/local/sw/SPAdes-3.14.0/bin/* /usr/local/bin

# FastQC
RUN cd /tmp/build && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod 775  FastQC/fastqc && \
    mv FastQC /usr/local/sw && \
    ln -s /usr/local/sw/FastQC/fastqc /usr/local/bin

# multiqc
RUN cd /tmp/build && \
    git clone https://github.com/ewels/MultiQC.git && \
    cd MultiQC && git checkout v1.8 && \
    python setup.py install

# add
# https://github.com/broadinstitute/picard/releases/download/2.22.0/picard.jar
# samtools sort -@ 2 -m 1G Aligned.out.bam -o aligned_sort.bam &> sam_sort.log &

#RUN git clone https://github.com/

#git clone https://github.com/
#cd workflows
#conda env create -f env/rnaseq.yml

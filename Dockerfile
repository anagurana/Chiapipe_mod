################################
## Draft Dockerfile to ChIA-PIPE
#
# Step 1. Install ChIA-PET Utilities (CPU)
# Step 2. Install other dependencies
# Step 3. Copy over shell scripts
#
# 28 Feb 2019 
##################################

# Set the base image to CentOS6
FROM centos

# File Author / Maintainer
MAINTAINER Daniel Capurso

# Install base packages: git, python, wget, unzip, R

RUN yum -y update && yum -y install \
    curl \
    git \
    tar \
    unzip \
    make \
    automake \
    gcc \
    gcc-c++ \
    kernel-devel \
    wget \
    bzip2 \
    zlib-devel.x86_64 \
    bc \
    #fot samtools
    ncurses-devel \
    bzip2-devel \
    #for R
    readline-devel \
    gcc-gfortran


RUN git clone https://github.com/anagurana/Chiapipe_mod.git

### Install ChIA-PET Utilities (CPU)

### Download and compile zlib (v 1.2.8)
#RUN curl -O https://www.zlib.net/fossils/zlib-1.2.8.tar.gz && tar -xzvf zlib-1.2.8.tar.gz
#RUN cd zlib-1.2.8 && ./configure && make
#
### Download and compile Chiasig
#RUN wget http://folk.uio.no/jonaspau/chiasig/ChiaSigCPPv093.zip && unzip ChiaSigCPPv093.zip
#RUN cd ChiaSigCPPv093 && make

### LOCAL INSTALL CHIA PIPE DEPENDENCIES

## Insert install directory at front of PATH
ARG install_dir=/Chiapipe_mod/dependencies
ENV PATH=${install_dir}:$PATH
WORKDIR ${install_dir}

## Install Anaconda2
RUN wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh &&\
    bash Anaconda2-5.2.0-Linux-x86_64.sh -b -p ${install_dir}/anaconda2 &&\
    ln -s anaconda2/bin/python python && ln -s anaconda2/bin/conda conda

## Install pysam, numpy, regex (Python packages), zMACS peak caller
RUN echo "Installing pysam..." && ./conda install -c anaconda pysam &&\
    echo "Installing regex..." && ./conda install -c anaconda regex &&\
    echo "Installing macs2..." && ./conda install -c bioconda macs2 && \
    ln -s anaconda2/bin/macs2 macs2

# Install biopython v 1.76
RUN wget http://biopython.org/DIST/biopython-1.76.tar.gz && \
    tar -xzvf biopython-1.76.tar.gz && \
    cd biopython-1.76 && python setup.py install

## Install pigz
RUN wget https://zlib.net/pigz/pigz.tar.gz && tar -xzvf pigz.tar.gz && \
    cd $(ls | grep "pigz-") && make && \
    cp pigz unpigz ../ && rm -r pigz

## Install java/1.8
RUN wget -c --header "Cookie: oraclelicense=accept-securebackup-cookie" \
    http://download.oracle.com/otn-pub/java/jdk/\
8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz && \
    tar -xzvf jdk-8u131-linux-x64.tar.gz && \
    cp -r jdk1.8.0_131/jre/ . && yes | rm -r jdk1.8.0_131 && \
    ln -s jre/bin/java java

## Install perl/5.26.0
RUN wget http://www.cpan.org/src/5.0/perl-5.26.0.tar.gz && \
    tar -xzvf perl-5.26.0.tar.gz && \
    cd perl-5.26.0 && ./Configure -des -Dprefix=${install_dir} && \
    make && make install && cd ../ && \
    yes | rm -r perl-5.26.0

## Install bedtools/2.26.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/\
bedtools-2.26.0.tar.gz && tar -xzvf bedtools-2.26.0.tar.gz && \
    cd bedtools2 && make && \
    cd bin && cp * ../../ && \
    cd ../../ && rm -r bedtools2

## Install samtools/1.5
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.5/\
samtools-1.5.tar.bz2 && tar -xvjf samtools-1.5.tar.bz2 && \
    cd samtools-1.5 && ./configure --disable-lzma && make && \
    cp samtools ../ && cd ../ && rm -r samtools-1.5

## Install R/3.2.1
RUN wget http://lib.stat.cmu.edu/R/CRAN/src/base/R-3/R-3.2.1.tar.gz && \
    tar -xzvf R-3.2.1.tar.gz && cd R-3.2.1 && \
    ./configure --prefix=${install_dir} --with-x=no && make && \
    ln -s R-3.2.1/bin/R R

# ChIA-PIPE excecution
WORKDIR /Chiapipe_mod
COPY fastq fastq
COPY reference reference
COPY config_file.sh config_file.sh


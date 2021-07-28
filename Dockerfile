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

## Dependencies
# zlib v 1.2.8
# ChiaSigScaled v 19.4.17n
# bwa       (Compiled during CPU make ??)
# BWT-SW    (Compiled during CPU make ??)

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
    bzip2

RUN git clone https://github.com/anagurana/Chiapipe_mod.git
WORKDIR /Chiapipe_mod

# Create directory for the project
# RUN mkdir chia-pipe

### Install ChIA-PET Utilities (CPU)

## Make outer working directory
#WORKDIR chia-pipe/util/cpu-dir
#
### Download and compile zlib (v 1.2.8)
#RUN curl -O https://www.zlib.net/fossils/zlib-1.2.8.tar.gz && tar -xzvf zlib-1.2.8.tar.gz
#RUN cd zlib-1.2.8 && ./configure && make
#
### Download and compile Chiasig
#RUN wget http://folk.uio.no/jonaspau/chiasig/ChiaSigCPPv093.zip && unzip ChiaSigCPPv093.zip
#RUN cd ChiaSigCPPv093 && make
#
### Download and compile bwa
#RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 && tar -xvf bwa-0.7.17.tar.bz2
#RUN cd bwa-0.7.17 && make

### LOCAL INSTALL CHIA PIPE DEPENDENCIES

## Insert install directory at front of PATH
RUN install_dir=/Chiapipe_mod/dependencies
WORKDIR ${install_dir}
RUN export PATH="${install_dir}:${PATH}"

## Install Anaconda2
RUN wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh
RUN bash Anaconda2-5.2.0-Linux-x86_64.sh -b -p ${install_dir}/anaconda2
RUN ln -s /anaconda2/bin/python python && ln -s /anaconda2/bin/conda conda

## Install pysam, biopython, numpy, regex (Python packages)
RUN echo "Installing pysam..." && ./conda install -c anaconda pysam
RUN echo "Installing biopython..." && ./conda install -c anaconda biopython
RUN echo "Installing numpy..." && ./conda install -c anaconda numpy
RUN echo "Installing regex..." && ./conda install -c anaconda regex

## Install MACS peak caller
RUN echo "Installing macs2..." && ./conda install -c bioconda macs2
RUN ln -s anaconda2/bin/macs2 macs2

WORKDIR /Chiapipe_mod/util/cpu-dir
RUN cd zlib-1.2.8 && ./configure && make

## Install pigz
RUN wget https://zlib.net/pigz/pigz.tar.gz && tar -xzvf pigz.tar.gz
RUN cd $(ls | grep "pigz-") && make && cp pigz unpigz ../
RUN rm -r pigz

## Install java/1.8
RUN wget -c --header "Cookie: oraclelicense=accept-securebackup-cookie" \
http://download.oracle.com/otn-pub/java/jdk/\
8u131-b11/d54c1d3a095b4ff2b6607d096fa80163/jdk-8u131-linux-x64.tar.gz && tar -xzvf jdk-8u131-linux-x64.tar.gz
RUN cp -r jdk1.8.0_131/jre/ . && yes | rm -r jdk1.8.0_131
RUN ln -s jre/bin/java java

## Install perl/5.26.0
RUN wget http://www.cpan.org/src/5.0/perl-5.26.0.tar.gz && tar -xzvf perl-5.26.0.tar.gz
RUN cd perl-5.26.0 && ./Configure -des -Dprefix=${install_dir}
RUN cd perl-5.26.0 && make && make test && make install
RUN yes | rm -r perl-5.26.0

## Install bedtools/2.26.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/\
bedtools-2.26.0.tar.gz && tar -xzvf bedtools-2.26.0.tar.gz
RUN cd bedtools2 && make && cd bin && cp * ../../ && cd ../../ && rm -r bedtools2

## Install samtools/1.5
RUN wget https://sourceforge.net/projects/samtools/files/samtools/1.5/\
samtools-1.5.tar.bz2 && tar -xvjf samtools-1.5.tar.bz2
RUN cd samtools-1.5 && ./configure --disable-lzma
RUN cd samtools-1.5 && make && cp samtools ../
RUN rm -r samtools-1.5

## Install R/3.2.1
RUN wget http://lib.stat.cmu.edu/R/CRAN/src/base/R-3/R-3.2.1.tar.gz && tar -xzvf R-3.2.1.tar.gz
RUN cd R-3.2.1 && ./configure --prefix=${install_dir} --with-x=no && make
RUN ln -s R-3.2.1/bin/R R

# Exit from dependencies directory
WORKDIR ../

 # COPY local_install_chia_pipe_dependencies.sh local_install_chia_pipe_dependencies.sh
 # RUN mkdir dependencies
 # RUN bash local_install_chia_pipe_dependencies.sh -i dependencies



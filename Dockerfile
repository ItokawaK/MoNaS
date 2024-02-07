FROM ubuntu:20.04 AS build
LABEL version="2.0"
LABEL maintainer="https://github.com/ItokawaK/MoNaS"

ARG TIMEZONE
RUN ln -snf /usr/share/zoneinfo/${TIMEZONE} /etc/localtime && echo ${TIMEZONE} > /etc/timezone

# For internal tests
RUN apt update && apt install -y ca-certificates && \
   sed -i 's/http:\/\/archive.ubuntu.com\/ubuntu\//https:\/\/linux.yz.yamagata-u.ac.jp\/ubuntu/' /etc/apt/sources.list

RUN apt-get update -y && \
   apt-get install -y --no-install-recommends \
   build-essential \
   wget \
   zlib1g-dev \
   zip unzip && \
   apt-get clean && \
   rm -rf /var/lib/apt/lists/*

RUN mkdir /opt/static

# Install samtools (The MIT License)
RUN cd /opt && \
  wget  https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 -O samtools.tar.bz2 && \
  mkdir samtools && \
  tar xjf samtools.tar.bz2 -C samtools --strip-components 1 && \
  (cd samtools && ./configure --without-curses --disable-bz2 --disable-lzma && make -j 4) && \
  cp samtools/samtools static && \
  rm samtools.tar.bz2

# Install bcftools (The MIT License)
RUN cd /opt && \
  wget  https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 -O bcftools.tar.bz2 && \
  mkdir bcftools && \
  tar xjf bcftools.tar.bz2 -C bcftools --strip-components 1 && \
  (cd bcftools && ./configure --disable-bz2 --disable-lzma && make -j 4) && \
  cp  bcftools/bcftools static && \
  rm  bcftools.tar.bz2

# Install bwa (GPL)
RUN cd /opt && \
  wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 -O bwa.tar.bz2 && \
  mkdir bwa && \
  tar xjf bwa.tar.bz2 -C bwa --strip-components 1 && \
  (cd bwa && make -j 4) && \
  cp bwa/bwa static && \
  rm bwa.tar.bz2

# Install bedtools (The MIT License)
RUN cd /opt/static && \
  wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static -O bedtools && \
  chmod +x bedtools

# Install freebayes (The MIT License)
RUN cd /opt/static && \
  wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz -O freebayes.gz && \
  gunzip freebayes.gz && \
  chmod +x freebayes

# Install muscle (Public domain)
RUN cd /opt/static && \
  wget https://drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux32.tar.gz && \
  tar xzvf muscle3.8.31_i86linux32.tar.gz && \
  mv muscle3.8.31_i86linux32 muscle && \
  rm muscle3.8.31_i86linux32.tar.gz

FROM ubuntu:20.04 AS release
ARG TIMEZONE
RUN ln -snf /usr/share/zoneinfo/${TIMEZONE} /etc/localtime && echo ${TIMEZONE} > /etc/timezone

RUN apt update && apt install -y ca-certificates && \
   sed -i 's/http:\/\/archive.ubuntu.com\/ubuntu\//https:\/\/linux.yz.yamagata-u.ac.jp\/ubuntu/' /etc/apt/sources.list

RUN apt-get update -y && \
   apt-get install -y --no-install-recommends python3.9 && \
   ln -s /usr/bin/python3.9 /usr/bin/python

RUN apt-get update -y && \
   apt-get install -y --no-install-recommends \
   python3-pip && \
   apt-get clean && \
   rm -rf /var/lib/apt/lists/* && \
   pip3 install biopython

COPY --from=build /opt/static/* /usr/local/bin
RUN mkdir /opt/MoNaS
COPY . /opt/MoNaS

RUN cd /opt/MoNaS && \
  python setup.py install && \
  rm -r /opt/*

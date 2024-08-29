FROM ubuntu:22.04

LABEL maintainer="daniel.podlesny@gmail.com"
LABEL version="1.2024.08"
LABEL description="samestr docker image"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev kraken2 seqtk
RUN apt clean

RUN mkdir -p /opt/software && cd /opt/software

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/software/miniconda3 && \
rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN git clone https://github.com/danielpodlesny/samestr.git && \
cd samestr && \  
  /opt/software/miniconda3/bin/conda install -y python=3.9 && \ 
  /opt/software/miniconda3/bin/conda env update -n base --file environment.yml && \
  /opt/software/miniconda3/bin/pip install . 

ENV PATH /opt/software/miniconda3/bin:$PATH

  

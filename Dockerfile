FROM condaforge/miniforge3:24.7.1-0

LABEL maintainer="daniel.podlesny@gmail.com"
LABEL version="1.2025.10"
LABEL description="samestr docker image"


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev kraken2 seqtk
RUN apt clean

ADD LICENSE MANIFEST.in README.md environment.yml pyproject.toml requirements.txt setup.py /opt/software/samestr/
ADD samestr /opt/software/samestr/samestr/

RUN cd /opt/software/samestr && \
  conda install -y --override-channels -c conda-forge -c bioconda python=3.9 && \
  conda env update -n base --file environment.yml
  
RUN cd /opt/software/samestr && pip install .

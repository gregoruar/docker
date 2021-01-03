FROM continuumio/miniconda3

RUN apt update && apt install -y git && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /tools && 
    git clone --recursive -C /tools https://github.com/nh13/DWGSIM &&
    cd /tools/DWGSIM &&
    git submodule init &&
    git submodule update &&
    make
   
ENV PATH /tools/DWGSIM:$PATH

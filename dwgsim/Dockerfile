FROM continuumio/miniconda3

RUN apt update && \ 
    apt install -y git make gcc libz-dev libncurses-dev && \
    apt-get clean && \ 
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /tools && \ 
    git clone --recursive https://github.com/nh13/DWGSIM /tools/DWGSIM && \
    cd /tools/DWGSIM && \
    git submodule init && \
    git submodule update && \
    make
   
ENV PATH /tools/DWGSIM:$PATH

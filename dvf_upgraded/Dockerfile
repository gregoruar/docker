# taken from mult1fractal/bx_DeepVirfinder
FROM continuumio/miniconda3
ENV VERSION 'no version available'
ENV TOOLNAME DeepVirFinder

RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --add channels default

RUN conda install python=3.6 numpy tensorflow=1.13 theano keras scikit-learn Biopython mkl-service fire && conda clean -a
RUN git clone https://github.com/jessieren/$TOOLNAME
RUN git clone https://github.com/gregoruar/docker_recipes.git
RUN chmod +x /$TOOLNAME/dvf.py /docker_recipes/dvf_upgraded/*.py
ENV PATH /$TOOLNAME:/docker_recipes/dvf_upgraded:$PATH 

FROM bioconductor/bioconductor_docker

WORKDIR /home/rstudio

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org')); BiocManager::install(ask=FALSE)"

RUN Rscript -e "options(repos = c(CRAN = 'https://cran.r-project.org'))"

RUN Rscript -e 'BiocManager::install("vitkl/ParetoTI", dependencies = c("Depends", "Imports", "LinkingTo"))'

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV SEGTOOLS_VERSION="1.2.4"

# libkrb5-dev needed for bigWigToBedGraph generation during genomedata creation to work
RUN apt-get update && apt-get install -y \
    bzip2 \
    libkrb5-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-4.7.12.1-Linux-x86_64.sh	 -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

RUN conda create -n reticulate_PCHA python=3.7.3 pip    

RUN Rscript -e "BiocManager::install('cowplot')"
RUN Rscript -e "BiocManager::install('archetypes')"
RUN Rscript -e "BiocManager::install('circlize')"
RUN Rscript -e "BiocManager::install('clusterProfiler')"
RUN Rscript -e "BiocManager::install('org.Hs.eg.db')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"

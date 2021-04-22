FROM ubuntu:16.04

MAINTAINER matsumoto <matsumoto@gen-info.osaka-u.ac.jp>

ENV VERSION 0.1.1

RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial-cran35/" | tee -a /etc/apt/sources.list \
 && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
 && apt-get update \
 && apt-get -y install build-essential git r-base-core libssl-dev libcurl4-openssl-dev automake autoconf\
 && cd ~ \
 && git clone https://github.com/lh3/minimap2 \
 && git clone https://github.com/lh3/bwa \
 && git clone --recursive https://github.com/samtools/htslib \
 && git clone https://github.com/samtools/samtools \
 && git clone https://github.com/ymatsumoto/mlstverse \
 && git clone https://github.com/ymatsumoto/mlstverse.Mycobacterium.db \
 && cd ~/minimap2 && make && cp minimap2 /usr/bin/ \
 && cd ~/bwa && make && cp bwa /usr/bin/ \
 && cd ~/htslib && autoheader && autoconf && ./configure && make && make install \
 && cd ~/samtools && autoheader && autoconf && ./configure && make && make install \
 && cd ~ \
 && echo "install.packages(\"BiocManager\", repos=\"https://cloud.r-project.org\")" >> install_mlstverse.r \
 && echo "BiocManager::install(\"Rsamtools\")" >> install_mlstverse.r \
 && echo "install.packages(\"devtools\", repos=\"https://cloud.r-project.org\")" >> install_mlstverse.r \
 && echo "devtools::install(\"mlstverse\")" >> install_mlstverse.r \
 && echo "devtools::install(\"mlstverse.Mycobacterium.db\")" >> install_mlstverse.r \
 && Rscript install_mlstverse.r \
 && cd ~ && rm -rf minimap2 bwa htslib samtools mlstverse mlstverse.Mycobacterium.db install_mlstverse.r \
 && apt-get -y remove build-essential \
 && apt-get -y autoremove \
 && apt-get clean all \
 && groupadd -c -g 101 gen-info \
 && useradd -u 10466 nanopore-user 
ENTRYPOINT ["bash"]

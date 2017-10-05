FROM openjdk:8-jre

# Get and compile samtools and bowtie
RUN apt-get update && apt-get install -y build-essential zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev python libtbb2 \
    && wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2 \
    && tar xjvf samtools-1.4.1.tar.bz2 && rm samtools-1.4.1.tar.bz2 && cd samtools-1.4.1 \    
    && ./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.4.1 \
    && make all plugins-htslib \
    && cd / \ 
    && mv samtools-1.4.1 samtools \
    && wget "https://downloads.sourceforge.net/project/bowtie-bio/bowtie/1.0.0/bowtie-1.0.0-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie%2F1.0.0%2F&ts=1494353342&use_mirror=iweb" -O bowtie.zip \
    && unzip bowtie.zip && rm bowtie.zip \
    && mv bowtie-1.0.0 bowtie \
    && wget "https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.3/bowtie2-2.3.3-linux-x86_64.zip?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.3.3&ts=1504782363&use_mirror=kent" -O bowtie2.zip \
    && unzip bowtie2.zip && rm bowtie2.zip \
    && mv bowtie2-2.3.3 bowtie2 \
    && apt-get clean

ADD target/dist bicycle

ENV PATH "$PATH:/bowtie"
ENV PATH "$PATH:/bowtie2"
ENV PATH "$PATH:/samtools"
ENV PATH "$PATH:/bicycle/cmd"

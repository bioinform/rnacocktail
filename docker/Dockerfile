FROM ubuntu:16.04
		 

ENV SAMTOOLS_VERSION 1.2
ENV BEDTOOLS2_VERSION 2.24.0
ENV PYBEDTOOLS_VERSION 0.7.7
ENV PYSAM_VERSION 0.9.0
ENV HISAT2_VERSION 2.0.4
ENV STRINGTIE_VERSION 1.2.1
ENV OASES_VERSION 0.2.09
ENV VELVET_VERSION 1.2.10
ENV STAR_VERSION 2.5.2b
ENV PICARD_VERSION 2.2.2
ENV HTSLIB_VERSION 1.3
ENV GIREMI_VERSION 0.2.1
ENV FUSIONCATCHER_VERSION 0.99.7a

RUN apt-get update && \
    apt-get install -y --fix-missing build-essential zlib1g-dev unzip libncurses5-dev curl wget r-base r-base-dev python python-pip python-dev cmake libboost-all-dev  libxml2-dev libcurl4-gnutls-dev software-properties-common apt-transport-https default-jre default-jdk less vim

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
RUN apt-get update
RUN apt-get install -y --fix-missing r-base r-base-dev

RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})' > ~/.Rprofile
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2"); biocLite("tximport"); biocLite("readr");'

ADD https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2 
RUN cd /opt && tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && cd samtools-${SAMTOOLS_VERSION} && make && make install && cd /opt && rm -rf samtools*

ADD https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS2_VERSION}/bedtools-${BEDTOOLS2_VERSION}.tar.gz /opt/bedtools-${BEDTOOLS2_VERSION}.tar.gz
RUN cd /opt && tar -zxvf bedtools-${BEDTOOLS2_VERSION}.tar.gz && cd bedtools2 && make && make install && cd /opt && rm -rf bedtools*

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip -O /opt/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip  && cd /opt && unzip hisat2-${HISAT2_VERSION}-Linux_x86_64.zip && cp -p /opt/hisat2-${HISAT2_VERSION}/hisat2* /usr/local/bin && cd /opt && rm -rf hisat2*

ADD https://github.com/gpertea/stringtie/archive/v${STRINGTIE_VERSION}.tar.gz /opt/stringtie-${STRINGTIE_VERSION}.Linux_x86_64.tar.gz
RUN cd /opt && tar -zxvf stringtie-${STRINGTIE_VERSION}.Linux_x86_64.tar.gz && cd stringtie-${STRINGTIE_VERSION} && make && cp -p /opt/stringtie-${STRINGTIE_VERSION}/stringtie /usr/local/bin && cd /opt && rm -rf stringtie*

ADD https://github.com/COMBINE-lab/salmon/releases/download/v0.7.2/Salmon-0.7.2_linux_x86_64.tar.gz /opt/Salmon-0.7.2_linux_x86_64.tar.gz
RUN cd /opt && tar -zxvf Salmon-0.7.2_linux_x86_64.tar.gz && cp -p /opt/Salmon-0.7.2_linux_x86_64/bin/salmon /usr/local/bin && cp -p /opt/Salmon-0.7.2_linux_x86_64/lib/* /usr/local/lib && cd /opt && rm -rf Salmon*

ADD https://github.com/dzerbino/oases/archive/${OASES_VERSION}.tar.gz /opt/${OASES_VERSION}.tar.gz
RUN cd /opt && tar -zxvf ${OASES_VERSION}.tar.gz && rm -rf /opt/oases-${OASES_VERSION}/velvet /opt/${OASES_VERSION}.tar.gz

ADD https://www.ebi.ac.uk/~zerbino/velvet/velvet_${VELVET_VERSION}.tgz /opt/velvet_${VELVET_VERSION}.tgz
RUN cd /opt && tar -zxvf velvet_${VELVET_VERSION}.tgz && cd velvet_${VELVET_VERSION} && make OPENMP=1 && mv /opt/velvet_${VELVET_VERSION} /opt/oases-${OASES_VERSION}/velvet && cd /opt/oases-${OASES_VERSION} && make OPENMP=1 && cp -p /opt/oases-${OASES_VERSION}/oases /usr/local/bin && cp -p /opt/oases-${OASES_VERSION}/velvet/velvet* /usr/local/bin && rm -rf /opt/velvet_${VELVET_VERSION}.tgz
RUN rm -rf /opt/oases-${OASES_VERSION}/velvet/* && rm -rf /opt/oases-${OASES_VERSION}/velvet/.gitignore && rm -rf /opt/oases-${OASES_VERSION}/* && rm -rf /opt/oases*

RUN wget http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz -O /opt/subread-1.5.1-Linux-x86_64.tar.gz && cd /opt && tar -zxvf subread-1.5.1-Linux-x86_64.tar.gz  && cp -p /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin && cd /opt && rm -rf subread*

ADD https://github.com/GATB/gatb-core/archive/v1.1.0.tar.gz  /opt/gatb-core-1.1.0.tar.gz
RUN cd /opt && tar -zxvf gatb-core-1.1.0.tar.gz && cd gatb-core-1.1.0/gatb-core && mkdir build && cd build && cmake .. && make && make install && cd /opt && rm -rf gatb*

ADD http://www.atgc-montpellier.fr/download/sources/lordec/LoRDEC-0.6.tar.gz /opt/LoRDEC-0.6.tar.gz
RUN cd /opt && tar -zxvf LoRDEC-0.6.tar.gz && cd LoRDEC-0.6 && make && cp -p /opt/LoRDEC-0.6/lordec* /usr/local/bin && cd /opt && rm -rf LoRDEC*

ADD https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz /opt/STAR_${STAR_VERSION}.tar.gz
RUN cd /opt && tar -zxvf STAR_${STAR_VERSION}.tar.gz && cp -p /opt/STAR-${STAR_VERSION}/bin/Linux_x86_64_static/* /usr/local/bin && cd /opt && rm -rf STAR*

ADD https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip /opt/picard-tools-${PICARD_VERSION}.zip
RUN cd /opt && unzip picard-tools-${PICARD_VERSION}.zip  && cp -p /opt/picard-tools-${PICARD_VERSION}/* /usr/local/bin && cd /opt && rm -rf picard*

ADD https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 /opt/htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt && tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && cd htslib-${HTSLIB_VERSION} && ./configure && make && rm -rf /opt/htslib-${HTSLIB_VERSION}.tar.bz2


ADD https://github.com/zhqingit/giremi/archive/v${GIREMI_VERSION}.tar.gz /opt/giremi-${GIREMI_VERSION}.tar.gz
RUN cd /opt && tar -zxvf giremi-${GIREMI_VERSION}.tar.gz && cp -p giremi-${GIREMI_VERSION}/giremi* /usr/local/bin  && cd /opt && rm -rf giremi-*

RUN pip install --upgrade pip
RUN pip install pybedtools==${PYBEDTOOLS_VERSION} pysam==${PYSAM_VERSION} biopython==1.66 openpyxl==1.5.6 xlrd==0.6.1 numpy pandas

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download -O /opt/bowtie-1.1.2-linux-x86_64.zip
RUN cd /opt && unzip bowtie-1.1.2-linux-x86_64.zip && cp -p /opt/bowtie-1.1.2/bowtie* /usr/local/bin && cd /opt && rm -rf bowtie*

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download -O /opt/bowtie2-2.2.9-linux-x86_64.zip 
RUN cd /opt && unzip bowtie2-2.2.9-linux-x86_64.zip && cp -p /opt/bowtie2-2.2.9/bowtie2* /usr/local/bin && cd /opt && rm -rf bowtie2*

RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -O /opt/bwa-0.7.12.tar.bz2 
RUN cd /opt && tar xjf bwa-0.7.12.tar.bz2 && cd bwa-0.7.12 && make && cp -p /opt/bwa-0.7.12/bwa /usr/local/bin && cd /opt && rm -rf bwa*

ADD https://github.com/ndaniel/seqtk/archive/1.0-r82b.tar.gz /opt/seqtk-1.0-r82b.tar.gz
RUN cd /opt && tar -zxvf /opt/seqtk-1.0-r82b.tar.gz && cd seqtk-1.0-r82b && make  && cp -p /opt/seqtk-1.0-r82b/seqtk /usr/local/bin && cd /opt && rm -rf seqtk*

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat /usr/local/bin/blat
RUN chmod 755 /usr/local/bin/blat

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit /usr/local/bin/faToTwoBit
RUN chmod 755 /usr/local/bin/faToTwoBit

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v287/liftOver /usr/local/bin/liftOver
RUN chmod 755 /usr/local/bin/liftOver




ADD https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1/sratoolkit.2.8.1-ubuntu64.tar.gz /opt/sratoolkit.2.8.1-ubuntu64.tar.gz
RUN cd /opt && tar -zxvf sratoolkit.2.8.1-ubuntu64.tar.gz && cp -Rp /opt/sratoolkit.2.8.1-ubuntu64/bin/* /usr/local/bin/ && cd /opt && rm -rf sratoolkit*

ADD http://ftp.gnu.org/gnu/coreutils/coreutils-8.25.tar.xz /opt/coreutils-8.25.tar.xz 
RUN cd /opt && tar -xJf coreutils-8.25.tar.xz && cd coreutils-8.25 && ./configure FORCE_UNSAFE_CONFIGURE=1 && make && make install && cd /opt && rm -rf coreutils*

ADD http://pkgs.fedoraproject.org/repo/pkgs/pigz/pigz-2.3.1.tar.gz/e803f8bc0770c7a5e96dccb1d2dd2aab/pigz-2.3.1.tar.gz /opt/pigz-2.3.1.tar.gz
RUN cd /opt && tar -zxvf pigz-2.3.1.tar.gz && cd pigz-2.3.1 && make && cp -p /opt/pigz-2.3.1/pigz /usr/local/bin && cd /opt && rm -rf pigz*

RUN wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 -O /opt/samtools-0.1.19.tar.bz2 && cd /opt && tar -xjvf samtools-0.1.19.tar.bz2 && cd samtools-0.1.19 && make && rm -rf /opt/samtools-0.1.19.tar.bz2

ADD http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-02-15.tar.gz /opt/gmap-gsnap-2017-02-15.tar.gz
RUN cd /opt && tar -zxvf gmap-gsnap-2017-02-15.tar.gz && cd gmap-2017-02-15 && ./configure && make && make install && cd /opt && rm -rf gmap*

RUN wget https://sourceforge.net/projects/fusioncatcher/files/fusioncatcher_v${FUSIONCATCHER_VERSION}.zip -O /opt/fusioncatcher_v${FUSIONCATCHER_VERSION}.zip && cd /opt && unzip fusioncatcher_v${FUSIONCATCHER_VERSION}.zip  && cp -p /opt/fusioncatcher_v${FUSIONCATCHER_VERSION}/bin/sam2psl.py /usr/local/bin
ENV PATH $PATH:/opt/fusioncatcher_v${FUSIONCATCHER_VERSION}/bin/
RUN echo '[paths]' > /opt/fusioncatcher_v${FUSIONCATCHER_VERSION}/configuration.cfg
RUN echo 'samtools = /opt/samtools-0.1.19/' >> /opt/fusioncatcher_v${FUSIONCATCHER_VERSION}/configuration.cfg
RUN rm -rf /opt/fusioncatcher_v${FUSIONCATCHER_VERSION}.zip 

ADD https://github.com/bioinform/IDP/archive/master.zip /opt/IDP-master.zip
RUN rm -rf /usr/local/bin/IDP* && cd /opt && unzip IDP-master.zip && rm -rf /opt/IDP-master.zip

ADD https://www.healthcare.uiowa.edu/labs/au/IDP-fusion/files/IDP-fusion_1.1.1.tar.gz /opt/IDP-fusion_1.1.1.tar.gz
RUN cd /opt && tar -xzvf IDP-fusion_1.1.1.tar.gz && rm -rf /opt/IDP-fusion_1.1.1.tar.gz

RUN cd /opt
RUN cd /opt
RUN pip install https://github.com/bioinform/rnacocktail/archive/fix_examples.zip

VOLUME /work_dir

# RUN pip install https://github.com/bioinform/RNACocktail/archive/v0.1.tar.gz



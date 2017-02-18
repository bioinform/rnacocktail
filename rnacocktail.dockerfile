FROM ubuntu
		 

ENV SAMTOOLS_VERSION 1.2
ENV BEDTOOLS2_VERSION 2.26.0
ENV PYBEDTOOLS_VERSION 0.7.7
ENV PYSAM_VERSION 0.9.0
ENV HISAT2_VERSION 2.0.4
ENV STRINGTIE_VERSION 1.2.1
ENV OASES_VERSION 0.2.09
ENV VELVET_VERSION 1.2.10
ENV STAR_VERSION 2.5.1b 
ENV IDP_VERSION 0.1.9
ENV PICARD_VERSION 2.2.2
ENV HTSLIB_VERSION 1.3
ENV GIREMI_VERSION 0.2.1
ENV FUSIONCATCHER_VERSION 0.99.3e

RUN apt-get update && \
    apt-get install -y build-essential zlib1g-dev unzip libncurses5-dev curl wget r-base r-base-dev python python-pip python-dev libboost-all-dev cmake


# RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" > /etc/apt/sources.list.d/cran.list && \
# 	apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9 && \
#     add-apt-repository ppa:marutter/rdev && \
# 	apt-get update && \
# 	apt-get upgrade &&\
# 	apt-get install -y r-base=${R_BASE_VERSION}*

# RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | tee -a /etc/apt/sources.list
# RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
# RUN gpg -a --export E084DAB9 | sudo apt-key add -
# RUN apt-get update
# RUN apt-get install r-base r-base-dev

ADD https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 /opt/samtools-${SAMTOOLS_VERSION}.tar.bz2
RUN cd /opt && tar -xjvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && cd samtools-${SAMTOOLS_VERSION} && make && make install
RUN cd /opt && rm -rf samtools*

ADD https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS2_VERSION}/bedtools-${BEDTOOLS2_VERSION}.tar.gz /opt/bedtools-${BEDTOOLS2_VERSION}.tar.gz
RUN cd /opt && tar -zxvf bedtools-${BEDTOOLS2_VERSION}.tar.gz && cd bedtools2 && make && make install
RUN cd /opt && rm -rf bedtools*

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip -O /opt/hisat2-${HISAT2_VERSION}-Linux_x86_64.zip 
RUN cd /opt && unzip hisat2-${HISAT2_VERSION}-Linux_x86_64.zip
RUN cp -p /opt/hisat2-${HISAT2_VERSION}/hisat2-* /usr/local/bin
RUN cd /opt && rm -rf hisat2*

ADD http://ccb.jhu.edu/software/stringtie/dl/stringtie-${STRINGTIE_VERSION}.Linux_x86_64.tar.gz /opt/stringtie-${STRINGTIE_VERSION}.Linux_x86_64.tar.gz
RUN cd /opt && tar -zxvf stringtie-${STRINGTIE_VERSION}.Linux_x86_64.tar.gz 
RUN cp -p /opt/stringtie-${STRINGTIE_VERSION}.Linux_x86_64/stringtie /usr/local/bin
RUN cd /opt && rm -rf stringtie*

ADD https://github.com/COMBINE-lab/salmon/releases/download/v0.6.0/SalmonBeta-0.6.0_DebianSqueeze.tar.gz /opt/SalmonBeta-0.6.0_DebianSqueeze.tar.gz
RUN cd /opt && tar -zxvf SalmonBeta-0.6.0_DebianSqueeze.tar.gz
RUN cp -p /opt/SalmonBeta-0.6.1_DebianSqueeze/bin/salmon /usr/local/bin
RUN cd /opt && rm -rf Salmon*

ADD https://github.com/dzerbino/oases/archive/${OASES_VERSION}.tar.gz /opt/${OASES_VERSION}.tar.gz
RUN cd /opt && tar -zxvf ${OASES_VERSION}.tar.gz 
RUN rm -rf /opt/oases-${OASES_VERSION}/velvet

ADD https://www.ebi.ac.uk/~zerbino/velvet/velvet_${VELVET_VERSION}.tgz /opt/velvet_${VELVET_VERSION}.tgz
RUN cd /opt && tar -zxvf velvet_${VELVET_VERSION}.tgz && cd velvet_${VELVET_VERSION} && make OPENMP=1
RUN mv /opt/velvet_${VELVET_VERSION} /opt/oases-${OASES_VERSION}/velvet
RUN cd /opt/oases-${OASES_VERSION} && make OPENMP=1
RUN cp -p /opt/oases-${OASES_VERSION}/oases /usr/local/bin
RUN cp -p /opt/oases-${OASES_VERSION}/velvet/velvet* /usr/local/bin
RUN rm -rf /opt/oases-${OASES_VERSION}

RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cran.r-project.org"; options(repos=r)})' > ~/.Rprofile
RUN R -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2"); biocLite("tximport");'

ADD http://downloads.sourceforge.net/project/subread/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz /opt/subread-1.5.1-Linux-x86_64.tar.gz
RUN cd /opt && tar -zxvf subread-1.5.1-Linux-x86_64.tar.gz  
RUN cp -p /opt/subread-1.5.1-Linux-x86_64/bin/featureCounts /usr/local/bin
RUN cd /opt && rm -rf subread*

ADD https://github.com/GATB/gatb-core/archive/v1.1.0.tar.gz  /opt/gatb-core-1.1.0.tar.gz
RUN cd /opt && tar -zxvf gatb-core-1.1.0.tar.gz && cd gatb-core-1.1.0/gatb-core && mkdir build && cd build && cmake .. && make && make install
RUN cd /opt && rm -rf gatb*

ADD http://www.atgc-montpellier.fr/download/sources/lordec/LoRDEC-0.6.tar.gz /opt/LoRDEC-0.6.tar.gz
RUN cd /opt && tar -zxvf LoRDEC-0.6.tar.gz && cd LoRDEC-0.6 && make 
RUN cp -p /opt/LoRDEC-0.6/lordec* /usr/local/bin
RUN cd /opt && rm -rf LoRDEC*

ADD https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz /opt/STAR_${STAR_VERSION}.tar.gz
RUN cd /opt && tar -zxvf STAR_${STAR_VERSION}.tar.gz && cd STAR-${STAR_VERSION} && make STAR 
RUN cd /opt && rm -rf STAR*

ADD http://www.healthcare.uiowa.edu/labs/au/IDP/files/IDP_${IDP_VERSION}.tar.gz /opt/IDP_${IDP_VERSION}.tar.gz
RUN cd /opt && tar -zxvf IDP_${IDP_VERSION}.tar.gz && cp -R /opt/IDP_${IDP_VERSION}/ /usr/local/bin/IDP/
RUN cd /opt && rm -rf IDP*

ADD https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard-tools-${PICARD_VERSION}.zip /opt/picard-tools-${PICARD_VERSION}.zip
RUN cd /opt && unzip picard-tools-${PICARD_VERSION}.zip 
RUN cp -p /opt/picard-tools-${PICARD_VERSION}/* /usr/local/bin
RUN cd /opt && rm -rf picard*

ADD https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 /opt/htslib-${HTSLIB_VERSION}.tar.bz2
RUN cd /opt && tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && cd htslib-${HTSLIB_VERSION} && ./configure && make && make install
RUN cd /opt && rm -rf htslib*

ADD https://github.com/zhqingit/giremi/archive/v${GIREMI_VERSION}.tar.gz /opt/giremi-${GIREMI_VERSION}.tar.gz
RUN cd /opt && tar -zxvf giremi-${GIREMI_VERSION}.tar.gz && cp -p giremi-${GIREMI_VERSION}/giremi* /usr/local/bin 
RUN cd /opt && rm -rf giremi-*

ADD http://sf.net/projects/fusioncatcher/files/bootstrap.py /opt/bootstrap.py 
RUN python bootstrap.py -t --download

RUN pip install --upgrade pip
RUN pip install pybedtools==${PYBEDTOOLS_VERSION} pysam==${PYSAM_VERSION} biopython==1.66 openpyxl==1.5.6 xlrd==0.6.1 numpy pandas

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-linux-x86_64.zip/download -O /opt/bowtie-1.1.2-linux-x86_64.zip
RUN cd /opt && unzip bowtie-1.1.2-linux-x86_64.zip
RUN cp -p /opt/bowtie-1.1.2/bowtie* /usr/local/bin
RUN cd /opt && rm -rf bowtie*

RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip/download -O /opt/bowtie2-2.2.9-linux-x86_64.zip
RUN cd /opt && unzip bowtie2-2.2.9-linux-x86_64.zip
RUN cp -p /opt/bowtie2-2.2.9/bowtie2* /usr/local/bin
RUN cd /opt && rm -rf bowtie2*

RUN wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2/download -O /opt/bwa-0.7.12.tar.bz2
RUN cd /opt && tar xjf bwa-0.7.12.tar.bz2 && cd bwa-0.7.12 && make
RUN cp -p /opt/bwa-0.7.12/bwa /usr/local/bin
RUN cd /opt && rm -rf bwa*

ADD https://github.com/ndaniel/seqtk/archive/1.0-r82b.tar.gz /opt/seqtk-1.0-r82b.tar.gz
RUN cd /opt && tar -zxvf /opt/seqtk-1.0-r82b.tar.gz && cd seqtk-1.0-r82b && make 
RUN cp -p /opt/seqtk-1.0-r82b/seqtk /usr/local/bin
RUN cd /opt && rm -rf seqtk*

ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat /usr/local/bin/blat
ADD http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit /usr/local/bin/faToTwoBit

ADD https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1/sratoolkit.2.8.1-ubuntu64.tar.gz /opt/sratoolkit.2.8.1-ubuntu64.tar.gz
RUN cd /opt && tar -zxvf sratoolkit.2.8.1-ubuntu64.tar.gz
RUN cp -p /opt/sratoolkit.2.8.1-ubuntu64/bin/fastdump.2.8.1 /usr/local/bin
RUN cd /opt && rm -rf sratoolkit*

ADD http://ftp.gnu.org/gnu/coreutils/coreutils-8.25.tar.xz /opt/coreutils-8.25.tar.xz
RUN cd /opt && tar -xJf coreutils-8.25.tar.xz && cd coreutils-8.25 && ./configure FORCE_UNSAFE_CONFIGURE=1 && make && make install
RUN cd /opt && rm -rf coreutils*

ADD http://pkgs.fedoraproject.org/repo/pkgs/pigz/pigz-2.3.1.tar.gz/e803f8bc0770c7a5e96dccb1d2dd2aab/pigz-2.3.1.tar.gz /opt/pigz-2.3.1.tar.gz
RUN cd /opt && tar -zxvf pigz-2.3.1.tar.gz && cd pigz-2.3.1 && make
RUN cp -p /opt/pigz-2.3.1/pigz /usr/local/bin
RUN cd /opt && rm -rf pigz*

RUN wget http://sourceforge.net/projects/samtools/files/samtools/0.1.19/samtools-0.1.19.tar.bz2 -O /opt/samtools-0.1.19.tar.bz2
RUN cd /opt && tar -xjvf samtools-0.1.19.tar.bz2 && cd samtools-0.1.19 && make 

ADD https://github.com/ndaniel/fusioncatcher/archive/v${FUSIONCATCHER_VERSION}.zip /opt/fusioncatcher-${FUSIONCATCHER_VERSION}.zip
RUN cd /opt && unzip fusioncatcher-${FUSIONCATCHER_VERSION}.zip 
RUN cp -p /opt/fusioncatcher-${FUSIONCATCHER_VERSION}/fusioncatcher /usr/local/bin
RUN cp -p /opt/fusioncatcher-${FUSIONCATCHER_VERSION}/fusioncatcher-build /usr/local/bin
RUN cp -p /opt/fusioncatcher-${FUSIONCATCHER_VERSION}/fusioncatcher-batch /usr/local/bin

RUN pip install https://github.com/bioinform/RNACocktail/archive/v0.1.tar.gz

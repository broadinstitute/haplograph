FROM continuumio/miniconda3:24.1.2-0
ARG work_dir=/hla-la
WORKDIR ${work_dir}

RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        gawk \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev \
        python-dev \
        bsdmainutils \
        default-jre \
        default-jdk \
        unzip \
        python3-pycurl \
        python3-edlib \
        cpanminus \
        libboost-all-dev
        
# rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN cargo --help

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# BWA
RUN git clone https://github.com/lh3/bwa.git \
    && cd bwa \
    && make \
    && cp ./bwa /usr/bin

# HTSLIB
ARG htslib_version=1.19.1
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && mv /usr/local/bin/samtools /usr/bin/samtools \
    && samtools --help


# install picard
ARG picard_version=1.123
RUN wget https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard-tools-${picard_version}.zip \
    && unzip picard-tools-${picard_version}.zip \
    && mkdir -p /usr/local/bin/picard \
    && mv ${work_dir}/picard-tools-${picard_version}/* /usr/local/bin/picard/ \
    && chmod 0644 /usr/local/bin/picard/SamToFastq.jar \
    && rm -rf /usr/local/bin/picard-tools-${picard_version}/


# BCFTOOLS
ARG bcftools_version=1.21
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j 16 \
    && make install \
    && cd ${work_dir} \
    && bcftools --help

# install bamtools
ARG bamtools_version=2.5.1
RUN wget https://github.com/pezmaster31/bamtools/archive/v${bamtools_version}.zip \
    && unzip v${bamtools_version}.zip \
    && cd bamtools-${bamtools_version} \
    && mkdir build \
    && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=/usr/local/bin .. \
    && make \
    && make install


# Set environment variables for Bio::DB::HTS installation
# ENV HTSLIB_DIR=/usr/local/lib
# ENV HTSLIB_INC=/usr/local/include
# RUN cpanm --notest Text::LevenshteinXS List::MoreUtils
# RUN cpanm --notest IPC::Run XML::LibXML XML::Parser::PerlSAX XML::DOM XML::Twig
# RUN cpanm --notest Bio::SeqFeature::Lite

# # Install Bio::DB::HTS
# RUN cpanm --force Bio::DB::HTS
# RUN conda config --add channels bioconda \
#     && conda config --add channels conda-forge \
#     && conda config --set channel_priority strict \
#     && conda install hla-la --yes

# install HLA-LA
WORKDIR /usr/local/bin
ARG hla_la_version=1.0.4
RUN mkdir -p /usr/local/bin/HLA-LA/bin
RUN mkdir -p /usr/local/bin/HLA-LA/src
RUN mkdir -p /usr/local/bin/HLA-LA/obj
RUN mkdir -p /usr/local/bin/HLA-LA/temp
RUN mkdir -p /usr/local/bin/HLA-LA/working
RUN mkdir -p /usr/local/bin/HLA-LA/graphs
WORKDIR /usr/local/bin/
RUN wget https://github.com/DiltheyLab/HLA-LA/archive/refs/tags/v${hla_la_version}.zip
RUN unzip v${hla_la_version}.zip
RUN mv /usr/local/bin/HLA-LA-${hla_la_version}/* /usr/local/bin/HLA-LA/src/
RUN rm -rf /usr/local/bin/HLA-LA-${hla_la_version}
WORKDIR /usr/local/bin/HLA-LA/src
RUN sed -i 's@\$(BAMTOOLS_PATH)/lib64@\$(BAMTOOLS_PATH)/lib@' makefile
RUN make all BOOST_PATH=/usr/include/boost BAMTOOLS_PATH=/usr/local/bin

# modify paths.ini for hla*la
RUN sed -i 's@samtools_bin=@samtools_bin=/usr/local/bin/samtools@' paths.ini
RUN sed -i 's@bwa_bin=@bwa_bin=/usr/local/bin/bwa@' paths.ini
RUN sed -i 's@picard_sam2fastq_bin=.*@picard_sam2fastq_bin=/usr/local/bin/picard/SamToFastq.jar@' paths.ini

# some final cleanup and command
WORKDIR /usr/local/bin
ENV PATH="/usr/local/bin/HLA-LA/bin/:${PATH}"
ENV PATH="/usr/local/bin/HLA-LA/src/:${PATH}"
FROM continuumio/miniconda3:24.1.2-0
ARG work_dir=/truvari
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
        python3-edlib 
        

# GSUTIL
ARG gcloud_version=405.0.1
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

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
    && samtools --help

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


RUN pip install pipx
RUN pipx install virtualenv

ENV PATH=/root/.local/bin:${PATH}



# ENV PATH=${work_dir}/bin:${PATH}
# minigraph-cactus
ARG minigraphcactus_version=3.0.1
RUN wget https://github.com/ComparativeGenomicsToolkit/cactus/releases/download/v${minigraphcactus_version}/cactus-bin-v${minigraphcactus_version}.tar.gz \
    && tar -xzf cactus-bin-v${minigraphcactus_version}.tar.gz \
    && rm cactus-bin-v${minigraphcactus_version}.tar.gz \
    && cd cactus-bin-v${minigraphcactus_version} \
    && virtualenv -p python3 venv-cactus-v${minigraphcactus_version} \
    && printf "export PATH=$(pwd)/bin:\$PATH\nexport PYTHONPATH=$(pwd)/lib:\$PYTHONPATH\nexport LD_LIBRARY_PATH=$(pwd)/lib:\$LD_LIBRARY_PATH\n" >> venv-cactus-v${minigraphcactus_version}/bin/activate \
    && . venv-cactus-v${minigraphcactus_version}/bin/activate \
    && python3 -m pip install -U setuptools pip wheel \
    && python3 -m pip install -U . \
    && python3 -m pip install -U -r ./toil-requirement.txt \
    && cd bin && for i in wigToBigWig faToTwoBit bedToBigBed bigBedToBed axtChain pslPosTarget bedSort hgGcPercent mafToBigMaf hgLoadMafSummary hgLoadChain; do wget -q http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/${i}; chmod +x ${i}; done \
    && build-tools/downloadVCFWave \
    && cactus ./jobstore ./examples/evolverMammals.txt ./evolverMammals.hal

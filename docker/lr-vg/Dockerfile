# vg + htslib/samtools/bcftools + hifiasm for haplograph graph workflows.
# Base bumped off bullseye miniconda 24.1.2 (apt security mirror 404s).
FROM continuumio/miniconda3:25.1.1-2

ARG DEBIAN_FRONTEND=noninteractive
ARG work_dir=/vg
WORKDIR ${work_dir}

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses-dev \
        libffi-dev \
        libopenblas-dev \
        ca-certificates \
        gawk \
        tree \
        gnupg \
        openssh-client \
        time \
        curl \
        wget \
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
        python3-dev \
        xz-utils \
        tk-dev \
        bsdmainutils \
        default-jre-headless \
        unzip \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir -U crcmod edlib pycurl

# GSUTIL / gcloud
ARG gcloud_version=405.0.1
RUN wget -q https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && ./google-cloud-sdk/install.sh -q \
    && ln -s ${work_dir}/google-cloud-sdk/bin/gcloud /usr/local/bin/gcloud \
    && ln -s ${work_dir}/google-cloud-sdk/bin/gsutil /usr/local/bin/gsutil
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
ARG htslib_version=1.19.1
RUN wget -q https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make -j"$(nproc)" \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && ldconfig \
    && bgzip --help >/dev/null

# SAMTOOLS
ARG samtools_version=1.19.2
RUN wget -q https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make -j"$(nproc)" \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --version

# BCFTOOLS
ARG bcftools_version=1.21
RUN wget -q https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make -j"$(nproc)" \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --version

# VG (prebuilt static-ish binary from release)
ARG vg_version=1.77.0
RUN wget -q -O /usr/local/bin/vg \
        https://github.com/vgteam/vg/releases/download/v${vg_version}/vg \
    && chmod +x /usr/local/bin/vg \
    && vg version

# Hifiasm
ARG hifiasm_version=0.25.0
RUN wget -q https://github.com/chhylp123/hifiasm/archive/refs/tags/${hifiasm_version}.tar.gz \
    && tar xzf ${hifiasm_version}.tar.gz \
    && rm ${hifiasm_version}.tar.gz \
    && cd hifiasm-${hifiasm_version} \
    && make -j"$(nproc)" \
    && cp hifiasm /usr/local/bin/ \
    && cd ${work_dir} \
    && rm -rf hifiasm-${hifiasm_version} \
    && hifiasm --version

CMD ["/bin/bash"]

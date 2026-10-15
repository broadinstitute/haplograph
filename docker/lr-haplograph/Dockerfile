FROM rust:trixie
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt-get install -y --no-install-recommends git ssh build-essential make cmake zlib1g-dev vim wget unzip libclang-dev\
    && rm -rf /var/lib/apt/lists/*

RUN git config --global http.sslverify false
RUN git clone https://github.com/broadinstitute/haplograph.git \ 
    && cd haplograph \
    && cargo build --release
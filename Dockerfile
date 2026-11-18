# Flye Docker Image
# De novo assembler for single molecule sequencing reads
# https://github.com/fenderglass/Flye

FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build and runtime dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-pip \
    python3-dev \
    git \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /tmp/flye-build

# Copy source code
COPY . .

# Install Python dependencies and build Flye
RUN pip3 install --no-cache-dir setuptools && \
    python3 setup.py install && \
    rm -rf /tmp/flye-build

# Test installation
RUN flye --version

# Create working directory for data
WORKDIR /data

# Add labels
LABEL maintainer="Flye Developers"
LABEL description="Flye: de novo assembler for single molecule sequencing reads"
LABEL version="2.9.6"

# Set entrypoint
ENTRYPOINT ["flye"]
CMD ["--help"]

# Usage:
# Build the Docker image:
# docker build -t krona-kaiju .
#
# Run the container:
# docker run -it \
#     -v /path/to/data:/data \
#     krona-kaiju <input_R1.fastq.gz> <input_R2.fastq.gz> <output_prefix> <database_path>
#
# Example:
# docker run -it \
#     -v /Users/user/data:/data \
#     krona-kaiju sample_R1.fastq.gz sample_R2.fastq.gz sample_output /opt/kaijudb

# Base Image
FROM ubuntu:20.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update --allow-releaseinfo-change && \
    apt-get update && apt-get install -y \
        wget \
        tar \
        perl \
        gcc \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libwww-perl \
        libxml-simple-perl \
        libhtml-form-perl \
        libhtml-tree-perl \
        libio-socket-ssl-perl \
        libxml-libxml-perl \
        libxml-sax-expat-perl \
        libnet-http-perl \
        liblwp-protocol-https-perl \
        default-jre \
        nano && \
    apt-get clean

# Copy pre-downloaded Kaiju databases from the host
WORKDIR /opt/kaijudb
COPY kaiju_db_refseq_2024-08-14.tgz /opt/kaijudb/
COPY kaiju_db_viruses_2024-08-15.tgz /opt/kaijudb/

# Extract Kaiju databases
RUN tar -xvzf kaiju_db_refseq_2024-08-14.tgz && \
    tar -xvzf kaiju_db_viruses_2024-08-15.tgz && \
    rm *.tgz

# Install Kaiju
WORKDIR /opt
RUN wget https://github.com/bioinformatics-centre/kaiju/archive/refs/tags/v1.10.1.tar.gz && \
    tar -xvzf v1.10.1.tar.gz && \
    rm v1.10.1.tar.gz && \
    cd kaiju-1.10.1/src && \
    make && \
    ln -s /opt/kaiju-1.10.1/bin/kaiju* /usr/local/bin/

# Install KronaTools
WORKDIR /opt
RUN wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar && \
    tar -xvf KronaTools-2.8.1.tar && rm KronaTools-2.8.1.tar

# Add Kaiju and KronaTools to PATH
ENV PATH="/opt/KronaTools-2.8.1/scripts:/opt/kaiju-1.10.1/bin:$PATH"

# Set working directory
WORKDIR /data

# Copy the process script into the container
COPY process_krona-kaiju.sh /data/process_krona-kaiju.sh

# Make the script executable
RUN chmod +x /data/process_krona-kaiju.sh

# Set default command to run the script
CMD ["/data/process_krona-kaiju.sh"]


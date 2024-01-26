# Use Ubuntu as the base image
FROM ubuntu:latest

# Install dependencies and tools
RUN apt-get update && apt-get install -y \
    wget \
    ncbi-blast+

# Create a directory for the BLAST databases
RUN mkdir -p /blast/db

# Set the working directory
WORKDIR /blast/clustered-db

# Download the nr protein database
# Note: This is a large download and can take a significant amount of time
RUN wget -O nr_rep_seq.fasta.gz 'https://osf.io/gqxva/download' && \
    gunzip nr_rep_seq.fasta.gz

# Set environment variable for BLAST database
ENV BLASTDB /blast/clustered-db

# Default command on container start
CMD ["bash"]

# Base image for mmseqs2
FROM soedinglab/mmseqs2:latest

# Add multiple Debian mirrors for redundancy
RUN echo "deb http://deb.debian.org/debian bullseye main contrib non-free" > /etc/apt/sources.list && \
    echo "deb http://security.debian.org/debian-security bullseye-security main contrib non-free" >> /etc/apt/sources.list && \
    echo "deb http://mirror.csclub.uwaterloo.ca/debian/ bullseye main contrib non-free" >> /etc/apt/sources.list

# Install Python and Flask with retry logic
RUN apt-get update --allow-releaseinfo-change && \
    for i in {1..3}; do apt-get install -y python3 python3-pip && break || sleep 15; done

RUN pip3 install fastapi uvicorn python-multipart

# Copy the FastAPI app to the container
COPY app.py /usr/local/bin/app.py

# Set the working directory
WORKDIR /usr/local/bin

# Expose the port on which FastAPI will run
EXPOSE 8001

# Start the FastAPI server when the container starts
CMD ["python3", "app.py"]
# Base image for mmseqs2
FROM soedinglab/mmseqs2:latest

# Add the standard Debian repositories to ensure we can install all packages
RUN echo "deb http://deb.debian.org/debian bullseye main contrib non-free" > /etc/apt/sources.list
RUN echo "deb http://security.debian.org/debian-security bullseye-security main contrib non-free" >> /etc/apt/sources.list

# Install Python and Flask
RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install fastapi uvicorn

# Copy the FastAPI app to the container
COPY app.py /usr/local/bin/app.py

# Set the working directory
WORKDIR /usr/local/bin

# Expose the port on which FastAPI will run
EXPOSE 8001

# Start the FastAPI server when the container starts
CMD ["python3", "app.py"]

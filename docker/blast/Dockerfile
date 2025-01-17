FROM ncbi/blast

# Install Python and Flask
RUN apt-get update && apt-get install -y python3 python3-pip
RUN pip3 install fastapi uvicorn

# Add the Python script to the container
COPY app.py /usr/local/bin/app.py

# Set the working directory
WORKDIR /usr/local/bin

# Run the Python server script
CMD ["python3", "app.py"]

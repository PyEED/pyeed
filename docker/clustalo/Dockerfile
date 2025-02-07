FROM python:3.12

WORKDIR /app/api

# Install clustalo
RUN wget http://www.clustal.org/omega/clustalo-1.2.4-Ubuntu-x86_64 \
    && chmod +x clustalo-1.2.4-Ubuntu-x86_64 \
    && mv clustalo-1.2.4-Ubuntu-x86_64 /usr/local/bin/clustalo

# Install python dependencies
RUN pip install fastapi python-multipart uvicorn

COPY app.py .


CMD ["python", "app.py"]
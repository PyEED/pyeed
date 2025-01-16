sudo docker stop mmseqs
sudo docker remove mmseqs
sudo docker build --no-cache -t mmseqs_docker .
sudo docker run -d --name mmseqs -p 8001:8001 mmseqs_docker
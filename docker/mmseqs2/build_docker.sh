sudo docker stop mmseqs
sudo docker remove mmseqs

sudo docker build --no-cache -t mmseqs .
sudo docker run -d -p 8001:8001 --name mmseqs mmseqs
sudo docker stop clustalo
sudo docker remove clustalo

sudo docker build -t clustalo .
sudo docker run -d -p 6001:6001 clustalo

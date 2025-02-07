sudo docker stop clustalo
sudo docker remove clustalo

sudo docker build --no-cache -t clustalo .
sudo docker run -d -p 5001:5001 --name clustalo clustalo

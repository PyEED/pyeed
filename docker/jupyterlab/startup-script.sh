#!/bin/bash

# Download the latest requirements.txt from GitHub
echo "Downloading latest requirements.txt from GitHub..."
curl -o /tmp/requirements.txt -L https://raw.githubusercontent.com/PyEED/PyEED_JupyterLab/main/requirements.txt

# Update packages based on requirements.txt
echo "Checking for package updates..."
python /usr/local/bin/update_packages.py

# Start JupyterLab with no token/password (adjust as needed)
echo "Starting JupyterLab..."
exec start-notebook.sh --NotebookApp.token='' --NotebookApp.password=''

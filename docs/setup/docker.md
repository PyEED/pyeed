## Concept

The PyEED Docker Service combines the PyEED toolkit with JupyterLab, a web-based editor for writing and executing Jupyter Notebooks with analysis tools. In combination, the `pyeed` package can be used from inside JupyterLab, whereas the setup and configuration of the analysis tools are taken care of by the Docker Service(1).  
In this way, the Docker Service allows you to run JupyterLab on your local machine without having to install Python, Jupyter, and the necessary Python tools to work with your data.
{ .annotate }

1.  Docker is a way to run software in a container. This means that the software is isolated from the rest of your system. This simplifies the installation of computing environments since everything is preconfigured in the container. Docker is comparable to running a virtual machine, but instead of installing a whole operating system, it only installs the software you need.

To install Docker, follow the instructions for your operating system on the [Docker website](https://docs.docker.com/get-docker/).

## Initial Setup

1. **Install Docker**: Follow the instructions for your operating system on the [Docker website](https://docs.docker.com/get-docker/).

2. **[Download](https://github.com/PyEED/PyEED_JupyterLab/archive/refs/heads/main.zip
) the PyEED Docker Service**

4. **Start the Service** by running the following steps:
=== "Windows"

    1. Open the command line by pressing ++windows+r++ and type `powershell`.

    2. Navigate to the Downloads folder and unzip the downloaded file.

    3. Navigate to the unzipped folder by running the following command, adjust the path if necessary:
        ```powershell
        cd ~\Downloads\pyeed-main
        ```
    4. Start the Docker Service by running the following command:
        ```powershell
        docker compose up --build
        ```

=== "MacOS/Linux"

    1. Open the terminal
    2. Navigate to the Downloads folder and unzip the downloaded file.
        ```bash
        cd ~/Downloads
        unzip pyeed-main.zip
        ```
    3. Navigate to the unzipped folder by running the following command, adjust the path if necessary:
        ```bash
        cd ~/Downloads/pyeed-main
        ```
    4. Start the Docker Service by running the following command:
        ```bash
        docker compose up --build
        ```

## Start the PyEED Docker Service

After the initial setup, all containers belonging to the PyEED Docker Service are created and started. The service is now added to the `Containers` section in the Docker Desktop application. To start the service, click on the :material-play: button next to the `pyeed` container. To access the JupyterLab environment, click on the link `8888:8888` in the header of the container. This will open a new tab in your browser, showing the JupyterLab environment.

## Stopping the PyEED Docker Service

To stop the service, navigate to the `Containers` section in the Docker Desktop application and click on the :material-stop: button next to the `pyeed` container. Running containers are symbolized by a green container icon. You can close the browser window whenever you want. The container will keep running in the background unless you stop it in the Docker Desktop app.

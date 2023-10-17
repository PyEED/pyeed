# Setting up a Python environment and installing Docker

## Virtual environments

### Why should I manage environments?

Ideally, you should have one new virtual environment for every Python-based project you work on. The main purpose of this is to keep the dependencies of every project isolated from both the system and each other.

Find more information on why you should use virtual environments [here](https://blog.inedo.com/python/python-environment-management-best-practices/).

### Installing `miniconda`

Open this [link](https://docs.conda.io/projects/miniconda/en/latest/miniconda-install.html) and follow the instructions to install `miniconda` on your system.

To verify, that conda is installed properly, run the following command from a terminal:

```bash
conda list
```

### Creating a new environment

To create a new environment, run the following command from a terminal:

```bash
conda create --name pyeed_env python=3.10
```

In this case, the vrtual environment is named `pyeed_env` and uses Python 3.10.

### Activating the environment

To activate the environment, run the following command from a terminal:

```bash
conda activate pyeed_env
```

Now, every Python package you install will be installed in this environment and is isolated from the system and other environments.

### Installing pyEED

To install pyEED, run the following command from a terminal:

```bash
pip install git+https://github.com/PyEED/pyeed.git
```

## Docker

Docker is an open platform for developing, shipping, and running applications. Docker enables you to separate your applications from your infrastructure so you can deliver software quickly. With Docker, you can manage your infrastructure in the same ways you manage your applications. By taking advantage of Docker's methodologies for shipping, testing, and deploying code, you can significantly reduce the delay between writing code and running it in production.

For more information, visit the [Docker website](https://docs.docker.com/get-started/overview/).

### Installing Docker

Follow the instructions on the [Docker website](https://docs.docker.com/engine/install/) to install __Docker Desktop__ on your system.

### In case you use macOS

Currently, Docker has an issue with some mac version. Therefore, run the following command from a terminal:

```bash
export DOCKER_DEFAULT_PLATFORM=linux/amd64
```
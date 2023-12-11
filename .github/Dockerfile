# jupyter/scipy-notebook to build our image on top of it
FROM jupyter/scipy-notebook

LABEL Max Haeussler <max.haeussler@ibtb.uni-stuttgart.de>

USER root

# Install psycopg2 dependencies
RUN apt-get update \
    && apt-get -y install libpq-dev gcc python3-dev\
    && pip install psycopg2

# Install python libraries
RUN pip install --upgrade pip \
    && pip install git+https://github.com/PyEED/pyeed.git@alignments

# Let's change to  "$NB_USER" command so the image runs as a non root user by default
USER $NB_UID

#Let's define this parameter to install jupyter lab instead of the default juyter notebook command so we don't have to use it when running the container with the option -e
ENV JUPYTER_ENABLE_LAB=yes
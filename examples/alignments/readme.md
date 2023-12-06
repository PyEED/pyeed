# Create Docker image for Clustal Omega

The container is defined in the Dockerfile, which is used to build the image.

## Build the image

```bash
docker build -t clustal_omega .
```

## Run the container

The container can then be run within Python using `subprocess` module.

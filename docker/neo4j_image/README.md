## Neo4j Image Explanation

This is the Neo4j image for the Pyeed project. It is built using the Dockerfile in this directory.

The Dockerfile is based on the Neo4j Community Edition image.

When using the community edition, an issue can arise when wanting to perform a backup of the database.

The community edition does not support the `backup` command, so we need to use the enterprise edition for this. In order to perform the backup of the data folder the database needs to be stopped. This is not possible with the community edition. When do ing that the docker container stops as well so we can then obviolsy not perform the backup. As a workaround we can use the `neo4j-admin` command to perform the backup. But this only works if the docker container keeps running, this is possible if there is still an active command running in the container. This is why we use a cutsom image.

The custom image has a foreground process that keeps the container running. This foreground process is a `tail -f /dev/null` command. This command will keep the container running indefinitely.

When the container is started, the `neo4j start` command is executed to start the Neo4j server.


## Usage

A sample to start the container looks like this: (here the custom image is build at the name `my-neo4j`)
```bash

sudo docker run -d --name neo4j-niklas-example \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/password \
  -e NEO4J_ACCEPT_LICENSE_AGREEMENT="yes" \
  -e NEO4J_dbms_unmanaged__extension__classes="n10s.endpoint=/rdf" \
  -e NEO4J_dbms_security_procedures_unrestricted='["n10s.*", "apoc.*"]' \
  -e NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
  -e NEO4J_dbms_security_procedures_allowlist='["n10s.*", "apoc.*"]' \
  -p 1234:7687 \
  -p 4567:7474 \
  -v $(pwd)/data:/data \
  -v $(pwd)/import:/import \
  -v $(pwd)/plugins:/plugins \
  my-neo4j:latest
```


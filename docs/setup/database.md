# Neo4j Database

`pyeed` needs a database to store sequence data and analysis results. The database can be set up locally or on a hosted server.

## Local Setup

Docker can be used to set up a local Neo4j database. The following bash configures and starts a Neo4j database in a Docker container. The authentication and data directories can be adjusted.

```bash
docker run -it --name pyeed-neo4j \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/12345678 \
  -p 7474:7474 \
  -p 7687:7687 \
  -v $HOME/Documents/db/data:/data \
  -v $HOME/Documents/db/logs:/logs \
  -v $HOME/Documents/db/import:/var/lib/neo4j/import \
  -v $HOME/Documents/db/plugins:/plugins \
  -e NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="gds.*,apoc.*,gds.util.*" \
  -d neo4j:latest
```

This script configures the Neo4j database with enabled graph data science and apoc plugins.

## Hosted Setup

Neo4j offers a [free tier](https://neo4j.com/free-graph-database/) hosted database. If the free hosted database is inactive for a certain period, it will be deleted. To avoid data loss, it is recommended to back up the database regularly, use a local database, or upgrade to a paid plan.
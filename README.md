<div align="center">
<h1 align="center">pyeed</h1>
<h2 align="center">Python Enzyme Engineering Database</h2>
</div>


## About 📖
TBD

## Installation ⚙️

Install `pyeed` by running
```bash
pip install git+https://github.com/PyEED/pyeed.git
```

Setup Graph and Vector Databases via `docker-compose.yaml`:
```yaml
services:
  neo4j:
    image: neo4j:2025.10.1-community-bullseye
    container_name: neo4j
    restart: unless-stopped

    ports:
      - "7474:7474"   # Neo4j Browser / HTTP
      - "7687:7687"   # Bolt

    environment:
      NEO4J_AUTH: "neo4j/12345678"

      # memory tuning
      NEO4J_server_memory_heap_initial__size: 16g
      NEO4J_server_memory_heap_max__size: 16g
      NEO4J_server_memory_pagecache_size: 64g

      # plugins
      NEO4J_dbms_security_procedures_unrestricted: "apoc.*"
      NEO4J_PLUGINS: '["apoc", "graph-data-science"]'

    volumes:
      - /home/mha/dbs/neo4j/proteingraph/data:/data
      - /home/mha/dbs/neo4j/proteingraph/logs:/logs
      - /home/mha/dbs/neo4j/proteingraph/import:/import
      - /home/mha/dbs/neo4j/proteingraph/plugins:/plugins

    networks:
      - backend

  etcd:
    container_name: milvus-etcd
    image: quay.io/coreos/etcd:v3.5.18
    restart: unless-stopped
    environment:
      ETCD_AUTO_COMPACTION_MODE: revision
      ETCD_AUTO_COMPACTION_RETENTION: "1000"
      ETCD_QUOTA_BACKEND_BYTES: "4294967296"
      ETCD_SNAPSHOT_COUNT: "50000"
    command: >
      etcd
      -advertise-client-urls=http://etcd:2379
      -listen-client-urls=http://0.0.0.0:2379
      --data-dir=/etcd
    volumes:
      - /home/mha/dbs/milvus/etcd:/etcd
    healthcheck:
      test: ["CMD", "etcdctl", "endpoint", "health"]
      interval: 30s
      timeout: 20s
      retries: 3
    networks:
      - backend

  minio:
    container_name: milvus-minio
    image: minio/minio:RELEASE.2024-12-18T13-15-44Z
    restart: unless-stopped
    environment:
      MINIO_ROOT_USER: minioadmin
      MINIO_ROOT_PASSWORD: minioadmin
    command: ["server", "/minio_data", "--console-address", ":9001"]
    ports:
      - "9000:9000"   # S3 API
      - "9001:9001"   # MinIO console
    volumes:
      - /home/mha/dbs/milvus/minio:/minio_data
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9000/minio/health/live"]
      interval: 30s
      timeout: 20s
      retries: 3
    networks:
      - backend

  milvus:
    container_name: milvus-standalone
    image: milvusdb/milvus:v2.6.4
    restart: unless-stopped
    command: ["milvus", "run", "standalone"]
    security_opt:
      - seccomp:unconfined

    environment:
      ETCD_ENDPOINTS: etcd:2379
      MINIO_ADDRESS: minio:9000
      MQ_TYPE: woodpecker
      MINIO_ACCESS_KEY_ID: ${MINIO_ROOT_USER:-minioadmin}
      MINIO_SECRET_ACCESS_KEY: ${MINIO_ROOT_PASSWORD:-minioadmin}

    ports:
      - "19530:19530" # gRPC
      - "9091:9091"   # REST / metrics / healthz

    volumes:
      - /home/mha/dbs/milvus/data:/var/lib/milvus

    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:9091/healthz"]
      interval: 30s
      start_period: 90s
      timeout: 20s
      retries: 3

    depends_on:
      - etcd
      - minio
    networks:
      - backend

  attu:
    container_name: attu
    image: zilliz/attu:v2.6.1
    environment:
      MILVUS_URL: milvus-standalone:19530
    ports:
      - "8000:3000"
    depends_on:
      - "milvus"
    networks:
      - backend
      
networks:
  backend:
    driver: bridge

```

## Usage 🚀

...
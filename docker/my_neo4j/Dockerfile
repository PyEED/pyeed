FROM neo4j:latest

# Add a persistent foreground process
CMD ["sh", "-c", "/startup/docker-entrypoint.sh neo4j console && tail -f /dev/null"]
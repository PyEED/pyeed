from pyeed import Pyeed

def main(
    uri: str = "bolt://129.69.129.130:7687",
    user: str = "neo4j",
    password: str = "12345678",
) -> None:
    # Create a Pyeed object, automatically connecting to the database
    eedb = Pyeed(uri, user, password)

    eedb.db.generate_model_diagram(models_path="/home/nab/Niklas/pyeed/src/pyeed/model.py")

if __name__ == "__main__":
    main()
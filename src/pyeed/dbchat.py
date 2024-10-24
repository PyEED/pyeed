import openai
from neo4j.exceptions import CypherSyntaxError

from pyeed.dbconnect import DatabaseConnector


class DBChat:
    def __init__(self, db: DatabaseConnector) -> None:
        self.db = db

    def construct_cypher(
        self,
        question: str,
        openai_key: str,
        history=None,
    ):
        """Construct a Cypher query based on the provided question and history.
        Send the question to OpenAI's GPT-4 model to generate a Cypher query.

        Args:
            question (str): Question to ask the database.
            openai_key (str): OpenAI API key.
            history (optional): History of the conversation. Defaults to None.

        Returns:
            str: _description_
        """
        client = openai.OpenAI(api_key=openai_key)

        messages = [
            {"role": "system", "content": self.get_system_message()},
            {"role": "user", "content": question},
        ]

        # Used for Cypher healing flows
        if history:
            messages.extend(history)

        completions = client.chat.completions.create(
            messages=messages,  # type: ignore
            model="gpt-4",
            temperature=0.0,
            max_tokens=1000,
        )

        return completions.choices[0].message.content

    def get_system_message(self) -> str:
        return f"""
            Task: Generate Cypher queries to query a Neo4j graph database based on the provided schema definition.
            Instructions:
            Use only the provided relationship types and properties.
            Do not use any other relationship types or properties that are not provided.
            If you cannot generate a Cypher statement based on the provided schema, explain the reason to the user.
            Schema:
            {self.get_schema()}

            Note: Do not include any explanations or apologies in your responses.
            """

    def get_schema(self) -> str:
        """
        Returns the schema text for the PyEED database.
        """
        return f"""
            This is the schema representation of the Neo4j database.
            Node properties are the following:
            {self.db.node_properties}
            Relationship properties are the following:
            {self.db.relationship_properties}
            Relationship point from source to target nodes
            {self.db.relationships}
            Make sure to respect relationship types and directions.
            """

    def run(
        self,
        question: str,
        retry: bool,
        openai_key: str,
        history=None,
    ) -> list[dict]:
        """
        Query the database using natural language via OpenAI's GPT-4 model.

        Args:
            question (str): Question to ask the database.
            openai_key (str): OpenAI API key.
            retry (bool): Whether to retry once if the query if it fails.
        """
        # Construct Cypher statement
        cypher = self.construct_cypher(question, openai_key, history)
        try:
            return self.db.execute_read(cypher)
        # Self-healing flow
        except CypherSyntaxError as e:
            # If out of retries
            if not retry:
                return [{"error": str(e)}]
            # Self-healing Cypher flow by
            # providing specific error to GPT-4
            print("Retrying")
            return self.run(
                question,
                retry=False,
                openai_key=openai_key,
                history=[
                    {"role": "assistant", "content": cypher},
                    {
                        "role": "user",
                        "content": f"""This query returns an error: {str(e)} 
                        Give me a improved query that works without any explanations or apologies""",
                    },
                ],
            )

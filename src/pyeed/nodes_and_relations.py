from neomodel import (
    ArrayProperty,
    FloatProperty,
    IntegerProperty,
    StringProperty,
    StructuredNode,
    UniqueIdProperty,
    UniqueProperty,
)


class StrictStructuredNode(StructuredNode):
    """A StructuredNode subclass that raises an error if an invalid property is provided."""

    __abstract_node__ = True

    def __init__(self, *args, **kwargs):
        # Get the defined properties of the model
        allowed_properties = set(self.__class__._class_properties())

        # Check if any provided properties are not in the allowed set
        for key in kwargs:
            if key not in allowed_properties:
                raise AttributeError(
                    f"'{key}' is not a valid property for {self.__class__.__name__}"
                )

        super().__init__(*args, **kwargs)

    @classmethod
    def _class_properties(cls):
        """Retrieve all allowed properties (fields) defined on the class."""
        return {
            k
            for k, v in cls.__dict__.items()
            if isinstance(
                v,
                (
                    StringProperty,
                    IntegerProperty,
                    FloatProperty,
                    ArrayProperty,
                    UniqueIdProperty,
                ),
            )
        }

    def save(self, *args, **kwargs):
        """Validates the properties and then saves the node."""
        allowed_properties = self.__class__._class_properties()

        # Only validate properties defined in the model schema
        for field, prop in self.__dict__.items():
            if field not in allowed_properties:
                continue  # Skip non-class properties (like internal Neo4j fields)

            if prop is None or callable(prop):
                continue

            try:
                neo_type = getattr(self.__class__, field)
            except AttributeError:
                raise AttributeError(
                    f"'{self.__class__.__name__}' has no attribute '{field}'"
                )

            # Skip validation for UniqueIdProperty
            if isinstance(neo_type, UniqueIdProperty):
                continue

            # Validate StringProperty
            if isinstance(neo_type, StringProperty) and not isinstance(prop, str):
                raise TypeError(
                    f"Expected a string for '{field}', got {type(prop).__name__}"
                )

            # Validate IntegerProperty
            elif isinstance(neo_type, IntegerProperty) and not isinstance(prop, int):
                raise TypeError(
                    f"Expected an integer for '{field}', got {type(prop).__name__}"
                )

            # Validate FloatProperty
            elif isinstance(neo_type, FloatProperty) and not isinstance(prop, float):
                raise TypeError(
                    f"Expected a float for '{field}', got {type(prop).__name__}"
                )

            # Validate ArrayProperty
            elif isinstance(neo_type, ArrayProperty):
                if not isinstance(prop, list):
                    raise TypeError(
                        f"Expected a list for '{field}', got {type(prop).__name__}"
                    )

                # Validate list of integers, strings, or floats
                base_property = neo_type.base_property
                if isinstance(base_property, StringProperty):
                    if not all(isinstance(item, str) for item in prop):
                        raise TypeError(f"All items in '{field}' must be strings")
                elif isinstance(base_property, IntegerProperty):
                    if not all(isinstance(item, int) for item in prop):
                        raise TypeError(f"All items in '{field}' must be integers")
                elif isinstance(base_property, FloatProperty):
                    if not all(isinstance(item, float) for item in prop):
                        raise TypeError(f"All items in '{field}' must be floats")

        return super().save(*args, **kwargs)

    @classmethod
    def get_or_save(cls, **kwargs):
        """Attempts to save the node first, and if it already exists (due to unique constraint), retrieves it."""
        try:
            # Attempt to create and save a new node
            instance = cls(**kwargs)
            instance.save()
            return instance
        except UniqueProperty:
            # If a unique constraint error occurs, retrieve the existing node
            return cls.nodes.get(**kwargs)

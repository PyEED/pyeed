class TestCircular:
    def test_circular(self):
        """This test makes sure that there are no circular imports"""

        try:
            from pyeed import core
        except ImportError as e:
            assert False, f"Circular import detected:\n\n{e.msg}"

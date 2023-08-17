=======================
How to Write UI Testing
=======================

We use `pytest-qt <https://pytest-qt.readthedocs.io/en/latest/reference.html>`_
to emulate a user interacting with the GUI. This allows us to test the
GUI without having to manually click through the GUI.

See `tests/ui/` for examples.


Running Tests
=============

In order to run tests, in your terminal, navigate to the directory that
contains your tests files. There, you can simply use Pytest.

``pytest``

In order to run a specific test file you can just add the filename after
the pytest invocation:

``pytest file.py``


Other options include:

Multiple test files: ``pytest file1.py file2.py``

Single test: ``pytest -k "test_within_file"``

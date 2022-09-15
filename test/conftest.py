r"""
Fixtures for pytest
"""


# 3rd-party imports
import pytest

# standard imports
import os
import sys

this_module_path = sys.modules[__name__].__file__


@pytest.fixture(scope="module")
def data_server():
    r"""Object containing info and functionality for data files"""

    class _DataServe(object):

        _directory = os.path.join(os.path.dirname(this_module_path), "data")

        @property
        def directory(self):
            r"""Directory where to find the data es"""
            return self._directory

        def path_to(self, basename):
            r"""Absolute path to a data file. If it doesn't exist, try to find it in the remote repository"""
            file_path = os.path.join(self._directory, basename)
            for ext in ["", ".nxs.h5", "_event.nxs"]:
                if os.path.isfile(file_path + ext):
                    return file_path + ext
            raise IOError("File {0} not found in data directory {1}".format(basename, self._directory))

    return _DataServe()

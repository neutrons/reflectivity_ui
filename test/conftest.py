r"""
Fixtures for pytest
"""

from reflectivity_ui.interfaces.data_handling.filepath import RunNumbers
from reflectivity_ui.interfaces.data_handling.instrument import Instrument

# 3rd-party imports
import glob
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

        def get_file_paths(self, number):
            instrument = Instrument()

            run_numbers = RunNumbers(number)

            file_list = list()

            for run_number in run_numbers.numbers:
                search_string = instrument.file_search_template % run_number
                matches = glob.glob(search_string + ".nxs.h5")  # type: Optional[List[str]]
                if not matches:  # Look for old-style nexus file name
                    search_string = instrument.legacy_search_template % run_number
                    matches = glob.glob(search_string + "_event.nxs")
                if len(matches) >= 1:
                    file_list.append(matches[0])  # there should be only one match, since we query with one run number

            return file_list

    return _DataServe()

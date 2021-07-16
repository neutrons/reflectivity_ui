from __future__ import absolute_import, division, print_function, unicode_literals

# local imports
from reflectivity_ui.interfaces.data_handling import ApplicationConfiguration

# 3rd-party imports
import pytest

# standard imports
import os
import sys


class TestApplicationConfiguration(object):

    def test_init(self):
        # System's installation of mantid
        application_conf = ApplicationConfiguration()
        assert os.path.dirname(application_conf.mantid_path) in sys.path
        # Custom "installation" of mantid
        mantid_path = '/tmp/mantid41'
        if not os.path.isdir(mantid_path):
            os.makedirs('/tmp/mantid41')
        application_conf = ApplicationConfiguration(root_dir='/tmp')
        assert application_conf.mantid_path == '/tmp/mantid41'
        assert application_conf.mantid_version == '4.1.0'


if __name__ == '__main__':
    pytest.main([__file__])
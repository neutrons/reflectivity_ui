"""
    Data handling using Mantid
"""

import os


class ApplicationConfiguration(object):
    r"""Application-level configuration
    1. Selection of the mantid version and installation location:
        - Version is preferred in the following order, from most to least preferred: ['4.2.0', '4.1.0', '4.0.0']
        - Preferred installation location is under /opt, e.g. "/opt/mantid42". If no preferred installation is found,
          then a mantid installation is searched in the locations specified by `sys.path`.
          NOTE: the preferred installation is not tested for integrity. Only the existence of the location is tested,
          but not its contents.
    """
    # Polarization states
    POL_STATE = "PolarizerState"
    POL_VETO = "PolarizerVeto"
    ANA_STATE = "AnalyzerState"
    ANA_VETO = "AnalyzerVeto"

    _valid_mantid_versions = ['4.2.0', '4.1.0', '4.0.0']

    def __init__(self, root_dir='/opt'):
        self.mantid_path = None
        self.mantid_version = None
        # Check for installed Mantid versions under /opt
        for version in self._valid_mantid_versions:
            short = ''.join(version.split('.')[0:2])
            install_path = os.path.join(root_dir, 'mantid' + short)
            if os.path.isdir(install_path):  # check only for the existence of the directory, but not its contents
                self.mantid_path = install_path
                self.mantid_version = version
                return
        # Check for installed Mantid versions already in the PYTHONPATH
        try:
            import mantid
            self.mantid_path = os.path.dirname(mantid.__file__)
            self.mantid_version = mantid.__version__
        except ImportError:
            pass

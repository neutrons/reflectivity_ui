"""
    Data handling using Mantid
"""
class ApplicationConfiguration(object):
    """
        Application-level configuration
    """
    # Polarization states
    POL_STATE = "PolarizerState"
    POL_VETO = "PolarizerVeto"
    ANA_STATE = "AnalyzerState"
    ANA_VETO = "AnalyzerVeto"

    def __init__(self):
        self.mantid_path = '/opt/mantid42'

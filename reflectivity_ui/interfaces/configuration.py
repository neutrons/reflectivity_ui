from instrument import Instrument

class Configuration(object):
    """
        Hold reduction options
    """
    def __init__(self):
        self.instrument = Instrument()
        # Number of TOF bins
        self.tof_bins = 40
        # Bin type:
        #    0 = Constant bin width
        #    1 = Constant Q bin width
        #    2 = Cinstant 1/wavelength bin width
        self.tof_bin_type = 0
        ## Threshold under which we skip a cross-section, as fraction of the max count
        self.count_threshold = 0.01
    
    def to_q_settings(self):
        pass
    def from_q_settings(self):
        pass
    
class ApplicationConfiguration(object):
    def __init__(self):
        self.mantid_path = '.'
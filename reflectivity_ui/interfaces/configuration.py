"""
    Application configuration, including reduction options
"""
from __future__ import absolute_import, division, print_function
from .data_handling.instrument import Instrument

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
        #    2 = Constant 1/wavelength bin width
        self.tof_bin_type = 0

        # Threshold under which we skip a cross-section, as fraction of the max count
        self.count_threshold = 0.01
        self.tof_overwrite = None
        
        # Reduction parameters
        # Use region of interest specified in meta data
        self.use_roi = True

        # Update the specular peak range after finding the peak
        # within the ROI
        self.update_peak_range = False

        # Use background specified in the meta data, if available
        self.use_roi_bck = False
        self.use_tight_bck = False
        self.bck_offset = 0

        # Options to override the ROI, when available
        self.force_peak_roi = False
        self.forced_peak_roi = [0, 0]
        self.force_bck_roi = False
        self.forced_bck_roi = [0, 0]

    def to_q_settings(self):
        """ Save configuration to QSettings """
        pass
    def from_q_settings(self):
        """ Retrieve configuration from QSettings """
        pass

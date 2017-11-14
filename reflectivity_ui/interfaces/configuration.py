"""
    Application configuration, including reduction options
"""
from __future__ import absolute_import, division, print_function
import math

from .data_handling.instrument import Instrument

class Configuration(object):
    """
        Hold reduction options
    """
    def __init__(self, settings=None):
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
        self.bck_offset = 5

        # Options to override the range
        self.force_peak_roi = False
        self.forced_peak_roi = [120, 140]
        self.force_low_res_roi = False
        self.forced_low_res_roi = [120, 140]
        self.force_bck_roi = False
        self.forced_bck_roi = [10, 50]

        # Subtract background
        self.subtract_background = True
        # Overall scaling factor
        self.scaling_factor = 1.0
        # Cut first and last N points
        self.cut_first_n_points = 1
        self.cut_last_n_points = 1

        if settings is not None:
            self.from_q_settings(settings)

    def to_q_settings(self, settings):
        """
            Save configuration to QSettings
            :param settings QSettings: QSettings object
        """
        settings.setValue('use_roi', self.use_roi)
        settings.setValue('tof_bins', self.tof_bins)
        settings.setValue('tof_bin_type', self.tof_bin_type)
        settings.setValue('update_peak_range', self.update_peak_range)
        settings.setValue('use_roi_bck', self.use_roi_bck)
        settings.setValue('use_tight_bck', self.use_tight_bck)
        settings.setValue('bck_offset', self.bck_offset)

        settings.setValue('force_peak_roi', self.tof_bins)
        settings.setValue('forced_peak_roi', ','.join([str(x) for x in self.forced_peak_roi]))
        settings.setValue('force_low_res_roi', self.force_low_res_roi)
        settings.setValue('forced_low_res_roi', ','.join([str(x) for x in self.forced_low_res_roi]))
        settings.setValue('force_bck_roi', self.force_bck_roi)
        settings.setValue('forced_bck_roi', ','.join([str(x) for x in self.forced_bck_roi]))

        settings.setValue('subtract_background', self.subtract_background)
        settings.setValue('scaling_factor', self.scaling_factor)
        settings.setValue('cut_first_n_points', self.cut_first_n_points)
        settings.setValue('cut_last_n_points', self.cut_last_n_points)

    def from_q_settings(self, settings):
        """ Retrieve configuration from QSettings """
        self.use_roi = settings.value('use_roi', str(self.use_roi)).lower() == 'true'
        self.tof_bins = int(settings.value('tof_bins', self.tof_bins))
        self.tof_bin_type = int(settings.value('tof_bin_type', self.tof_bin_type))
        self.update_peak_range = settings.value('update_peak_range',
                                                str(self.update_peak_range)).lower() == 'true'
        self.use_roi_bck = settings.value('use_roi_bck',
                                          str(self.use_roi_bck)).lower() == 'true'
        self.use_tight_bck = settings.value('use_tight_bck',
                                            str(self.use_tight_bck)).lower() == 'true'
        self.bck_offset = int(settings.value('bck_offset', self.bck_offset))

        self.force_peak_roi = settings.value('force_peak_roi',
                                             str(self.force_peak_roi)).lower() == 'true'
        self.force_low_res_roi = settings.value('force_low_res_roi',
                                                str(self.force_low_res_roi)).lower() == 'true'
        self.force_bck_roi = settings.value('force_bck_roi',
                                            str(self.force_bck_roi)).lower() == 'true'

        default = ','.join([str(x) for x in self.forced_peak_roi])
        self.forced_peak_roi = [float(x) for x in settings.value('forced_peak_roi', default).split(',')]
        default = ','.join([str(x) for x in self.forced_low_res_roi])
        self.forced_low_res_roi = [float(x) for x in settings.value('forced_low_res_roi', default).split(',')]
        default = ','.join([str(x) for x in self.forced_bck_roi])
        self.forced_bck_roi = [float(x) for x in settings.value('forced_bck_roi', default).split(',')]

        self.subtract_background = settings.value('subtract_background', str(self.subtract_background)).lower() == 'true'
        self.scaling_factor = math.pow(10,float(settings.value('scaling_factor', self.scaling_factor)))
        self.cut_first_n_points = int(settings.value('cut_first_n_points', self.cut_first_n_points))
        self.cut_last_n_points = int(settings.value('cut_last_n_points', self.cut_last_n_points))


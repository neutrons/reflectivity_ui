"""
    Application configuration, including reduction options
"""
from __future__ import absolute_import, division, print_function

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
        self.set_direct_pixel = False
        self.direct_pixel_overwrite = 0
        self.set_direct_angle_offset = False
        self.direct_angle_offset_overwrite = 0
        self.use_dangle = False
        self.use_constant_q = False

        # Update the specular peak range after finding the peak
        # within the ROI
        self.update_peak_range = False

        # Use background specified in the meta data, if available
        self.use_roi_bck = False
        self.use_tight_bck = False
        self.bck_offset = 5

        # Options to override the range
        self.force_peak_roi = False
        self.peak_roi = [120, 140]
        self.force_low_res_roi = False
        self.low_res_roi = [120, 140]
        self.force_bck_roi = False
        self.bck_roi = [10, 50]

        # Subtract background
        self.subtract_background = True
        # Overall scaling factor
        self.scaling_factor = 1.0
        # Cut first and last N points
        self.cut_first_n_points = 1
        self.cut_last_n_points = 1

        # UI elements
        self.normalize_x_tof = False
        self.x_wl_map = False
        self.angle_map = False
        self.log_1d = True
        self.log_2d = True

        # Reduction options
        self.match_direct_beam = False
        self.normalization = None

        if settings is not None:
            self.from_q_settings(settings)

    @property
    def peak_position(self):
        return (self.peak_roi[1]+self.peak_roi[0])/2.0

    @peak_position.setter
    def peak_position(self, value):
        width = self.peak_width
        self.peak_roi[0] = value - width/2.0
        self.peak_roi[1] = value + width/2.0

    @property
    def peak_width(self):
        return self.peak_roi[1]-self.peak_roi[0]

    @peak_width.setter
    def peak_width(self, value):
        pos = self.peak_position
        self.peak_roi[0] = pos - value/2.0
        self.peak_roi[1] = pos + value/2.0

    @property
    def low_res_position(self):
        return (self.low_res_roi[1]+self.low_res_roi[0])/2.0

    @low_res_position.setter
    def low_res_position(self, value):
        width = self.low_res_width
        self.low_res_roi[0] = value - width/2.0
        self.low_res_roi[1] = value + width/2.0

    @property
    def low_res_width(self):
        return self.low_res_roi[1]-self.low_res_roi[0]

    @low_res_width.setter
    def low_res_width(self, value):
        pos = self.low_res_position
        self.low_res_roi[0] = pos - value/2.0
        self.low_res_roi[1] = pos + value/2.0

    @property
    def bck_position(self):
        return (self.bck_roi[1]+self.bck_roi[0])/2.0

    @bck_position.setter
    def bck_position(self, value):
        width = self.bck_width
        self.bck_roi[0] = value - width/2.0
        self.bck_roi[1] = value + width/2.0

    @property
    def bck_width(self):
        return self.bck_roi[1]-self.bck_roi[0]

    @bck_width.setter
    def bck_width(self, value):
        pos = self.bck_position
        self.bck_roi[0] = pos - value/2.0
        self.bck_roi[1] = pos + value/2.0

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
        settings.setValue('peak_roi', ','.join([str(x) for x in self.peak_roi]))
        settings.setValue('force_low_res_roi', self.force_low_res_roi)
        settings.setValue('low_res_roi', ','.join([str(x) for x in self.low_res_roi]))
        settings.setValue('force_bck_roi', self.force_bck_roi)
        settings.setValue('bck_roi', ','.join([str(x) for x in self.bck_roi]))

        settings.setValue('subtract_background', self.subtract_background)
        settings.setValue('scaling_factor', self.scaling_factor)
        settings.setValue('cut_first_n_points', self.cut_first_n_points)
        settings.setValue('cut_last_n_points', self.cut_last_n_points)

        settings.setValue('normalize_x_tof', self.normalize_x_tof)
        settings.setValue('x_wl_map', self.x_wl_map)
        settings.setValue('angle_map', self.angle_map)
        settings.setValue('log_1d', self.log_1d)
        settings.setValue('log_2d', self.log_2d)

        settings.setValue('use_constant_q', self.use_constant_q)
        settings.setValue('use_dangle', self.use_dangle)
        settings.setValue('set_direct_pixel', self.set_direct_pixel)
        settings.setValue('direct_pixel_overwrite', self.direct_pixel_overwrite)
        settings.setValue('set_direct_angle_offset', self.set_direct_angle_offset)
        settings.setValue('direct_angle_offset_overwrite', self.direct_angle_offset_overwrite)

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

        default = ','.join([str(x) for x in self.peak_roi])
        self.peak_roi = [float(x) for x in settings.value('peak_roi', default).split(',')]
        default = ','.join([str(x) for x in self.low_res_roi])
        self.low_res_roi = [float(x) for x in settings.value('low_res_roi', default).split(',')]
        default = ','.join([str(x) for x in self.bck_roi])
        self.bck_roi = [float(x) for x in settings.value('bck_roi', default).split(',')]

        self.subtract_background = settings.value('subtract_background', str(self.subtract_background)).lower() == 'true'
        self.scaling_factor = float(settings.value('scaling_factor', self.scaling_factor))
        self.cut_first_n_points = int(settings.value('cut_first_n_points', self.cut_first_n_points))
        self.cut_last_n_points = int(settings.value('cut_last_n_points', self.cut_last_n_points))

        self.normalize_x_tof = settings.value('normalize_x_tof', str(self.normalize_x_tof)).lower() == 'true'
        self.x_wl_map = settings.value('x_wl_map', str(self.x_wl_map)).lower() == 'true'
        self.angle_map = settings.value('angle_map', str(self.angle_map)).lower() == 'true'
        self.log_1d = settings.value('log_1d', str(self.log_1d)).lower() == 'true'
        self.log_2d = settings.value('log_2d', str(self.log_2d)).lower() == 'true'

        self.use_constant_q = settings.value('use_constant_q', str(self.use_constant_q)).lower() == 'true'
        self.use_dangle = settings.value('use_dangle', str(self.use_dangle)).lower() == 'true'
        self.set_direct_pixel = settings.value('set_direct_pixel', str(self.set_direct_pixel)).lower() == 'true'
        self.direct_pixel_overwrite = float(settings.value('direct_pixel_overwrite', self.direct_pixel_overwrite))
        self.set_direct_angle_offset = settings.value('set_direct_angle_offset', str(self.set_direct_angle_offset)).lower() == 'true'
        self.direct_angle_offset_overwrite = float(settings.value('direct_angle_offset_overwrite', self.direct_angle_offset_overwrite))



#pylint: disable=invalid-name, line-too-long, too-few-public-methods, too-many-instance-attributes, wrong-import-order, bare-except
"""
    Application configuration, including reduction options
"""
from __future__ import absolute_import, division, print_function
import sys
import logging
from .data_handling.instrument import Instrument

class Configuration(object):
    """
        Hold reduction options
    """
    # Choice of axes for off-specular binning
    QX_VS_QZ = 0
    KZI_VS_KZF = 1
    DELTA_KZ_VS_QZ = 3

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
        self.peak_position = 130
        self.peak_width = 20

        self.force_low_res_roi = False
        self.low_res_position = 130
        self.low_res_width = 20

        self.force_bck_roi = False
        self.bck_position = 30
        self.bck_width = 20

        # Subtract background
        self.subtract_background = True
        # Overall scaling factor
        self.scaling_factor = 1.0
        # Normalize to unity when stitching
        self.normalize_to_unity = True
        self.total_reflectivity_q_cutoff = 0.01

        # Cut first and last N points
        self.cut_first_n_points = 1
        self.cut_last_n_points = 1

        # UI elements
        self.normalize_x_tof = False
        self.x_wl_map = False
        self.angle_map = False
        self.log_1d = True
        self.log_2d = True

        # Off-specular options
        self.off_spec_x_axis = Configuration.DELTA_KZ_VS_QZ
        self.off_spec_slice = False
        self.off_spec_qz_list = []
        self.off_spec_err_weight = False
        self.off_spec_nxbins = 450
        self.off_spec_nybins = 200

        # Reduction options
        self.match_direct_beam = False
        self.normalization = None

        if settings is not None:
            try:
                self.from_q_settings(settings)
            except:
                logging.error("Could not process application settings\n  %s", sys.exc_value)

    @property
    def peak_roi(self):
        peak_min = int(round(self.peak_position - (self.peak_width - 1)/2.0))
        peak_max = int(round(self.peak_position + (self.peak_width - 1)/2.0))
        return [peak_min, peak_max]

    @peak_roi.setter
    def peak_roi(self, value):
        self.peak_position = (value[1] + value[0]) / 2.0
        self.peak_width = value[1] - value[0] + 1

    @property
    def low_res_roi(self):
        peak_min = int(round(self.low_res_position - (self.low_res_width - 1)/2.0))
        peak_max = int(round(self.low_res_position + (self.low_res_width - 1)/2.0))
        return [peak_min, peak_max]

    @low_res_roi.setter
    def low_res_roi(self, value):
        self.low_res_position = (value[1] + value[0]) / 2.0
        self.low_res_width = value[1] - value[0] + 1

    @property
    def bck_roi(self):
        peak_min = int(round(self.bck_position - (self.bck_width - 1)/2.0))
        peak_max = int(round(self.bck_position + (self.bck_width - 1)/2.0))
        return [peak_min, peak_max]

    @bck_roi.setter
    def bck_roi(self, value):
        self.bck_position = (value[1] + value[0]) / 2.0
        self.bck_width = value[1] - value[0] + 1

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

        # Off-specular options
        settings.setValue('off_spec_x_axis', self.off_spec_x_axis)
        settings.setValue('off_spec_slice', self.off_spec_slice)
        settings.setValue('off_spec_qz_list', ','.join([str(x) for x in self.off_spec_qz_list]))
        settings.setValue('off_spec_err_weight', self.off_spec_err_weight)
        settings.setValue('off_spec_nxbins', self.off_spec_nxbins)
        settings.setValue('off_spec_nybins', self.off_spec_nybins)

    def from_q_settings(self, settings):
        """ Retrieve configuration from QSettings """

        def _verify_true(parameter, default):
            """ Utility function to read a bool """
            _value = settings.value(parameter, str(default))
            return str(_value).lower() == 'true'

        self.use_roi = _verify_true('use_roi', self.use_roi)
        self.tof_bins = int(settings.value('tof_bins', self.tof_bins))
        self.tof_bin_type = int(settings.value('tof_bin_type', self.tof_bin_type))
        self.update_peak_range = _verify_true('update_peak_range', self.update_peak_range)
        self.use_roi_bck = _verify_true('use_roi_bck', self.use_roi_bck)
        self.use_tight_bck = _verify_true('use_tight_bck', self.use_tight_bck)
        self.bck_offset = int(settings.value('bck_offset', self.bck_offset))

        self.force_peak_roi = _verify_true('force_peak_roi', self.force_peak_roi)
        self.force_low_res_roi = _verify_true('force_low_res_roi', self.force_low_res_roi)
        self.force_bck_roi = _verify_true('force_bck_roi', self.force_bck_roi)

        default = ','.join([str(x) for x in self.peak_roi])
        self.peak_roi = [int(x) for x in settings.value('peak_roi', default).split(',')]
        default = ','.join([str(x) for x in self.low_res_roi])
        self.low_res_roi = [int(x) for x in settings.value('low_res_roi', default).split(',')]
        default = ','.join([str(x) for x in self.bck_roi])
        self.bck_roi = [int(x) for x in settings.value('bck_roi', default).split(',')]

        self.subtract_background = _verify_true('subtract_background', self.subtract_background)
        self.scaling_factor = float(settings.value('scaling_factor', self.scaling_factor))
        self.cut_first_n_points = int(settings.value('cut_first_n_points', self.cut_first_n_points))
        self.cut_last_n_points = int(settings.value('cut_last_n_points', self.cut_last_n_points))

        self.normalize_x_tof = _verify_true('normalize_x_tof', self.normalize_x_tof)
        self.x_wl_map = _verify_true('x_wl_map', self.x_wl_map)
        self.angle_map = _verify_true('angle_map', self.angle_map)
        self.log_1d = _verify_true('log_1d', self.log_1d)
        self.log_2d = _verify_true('log_2d', self.log_2d)

        self.use_constant_q = _verify_true('use_constant_q', self.use_constant_q)
        self.use_dangle = _verify_true('use_dangle', self.use_dangle)
        self.set_direct_pixel = _verify_true('set_direct_pixel', self.set_direct_pixel)
        self.direct_pixel_overwrite = float(settings.value('direct_pixel_overwrite', self.direct_pixel_overwrite))
        self.set_direct_angle_offset = _verify_true('set_direct_angle_offset', self.set_direct_angle_offset)
        self.direct_angle_offset_overwrite = float(settings.value('direct_angle_offset_overwrite', self.direct_angle_offset_overwrite))

        # Off-specular options
        self.off_spec_x_axis = int(settings.value('off_spec_x_axis', self.off_spec_x_axis))
        self.off_spec_slice = _verify_true('off_spec_slice', self.off_spec_slice)
        default = ','.join([str(x) for x in self.off_spec_qz_list])
        self.off_spec_qz_list = [float(x) for x in settings.value('off_spec_qz_list', default).split(',')]
        self.off_spec_err_weight = _verify_true('off_spec_err_weight', self.off_spec_err_weight)
        self.off_spec_nxbins = int(settings.value('off_spec_nxbins', self.off_spec_nxbins))
        self.off_spec_nybins = int(settings.value('off_spec_nybins', self.off_spec_nybins))

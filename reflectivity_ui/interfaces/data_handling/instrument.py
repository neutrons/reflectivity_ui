"""
    This instrument description contains information
    that is instrument-specific and abstracts out how we obtain
    information from the data file
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long
from __future__ import absolute_import, division, print_function
import sys
import math
import logging

# Import mantid according to the application configuration
from . import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *


# Constants
h = 6.626e-34  # m^2 kg s^-1
m = 1.675e-27  # kg

class Instrument(object):
    """
        Instrument class. Holds the data handling that is unique to a specific instrument.
    """
    n_x_pixel = 304
    n_y_pixel = 256
    huber_x_cut = 0
    peak_range_offset = 50
    tolerance = 0.05
    pixel_width = 0.0007
    instrument_name = "REFM"
    instrument_dir = "/SNS/REF_M"
    file_search_template = "/SNS/REF_M/*/data/REF_M_%s"

    def __init__(self):
        self.tof_range = [0,0]

    def get_tof_range(self, run_object):
        """
            Determine TOF range from the data
        """
        sample_detector_distance = run_object['SampleDetDis'].getStatistics().mean / 1000.0
        source_sample_distance = run_object['ModeratorSamDis'].getStatistics().mean / 1000.0
        source_detector_distance = source_sample_distance + sample_detector_distance

        wl = run_object.getProperty('LambdaRequest').value[0]
        chopper_speed = run_object.getProperty('SpeedRequest1').value[0]
        wl_offset = 0
        cst = source_detector_distance / h * m
        tof_min = cst * (wl + wl_offset * 60.0 / chopper_speed - 1.4 * 60.0 / chopper_speed) * 1e-4
        tof_max = cst * (wl + wl_offset * 60.0 / chopper_speed + 1.4 * 60.0 / chopper_speed) * 1e-4

        self.tof_range = [tof_min, tof_max]
        return [tof_min, tof_max]

    @classmethod
    def scattering_angle(cls, ws, peak_position=None):
        """
            Determine the scattering angle in degrees
        """
        dangle = ws.getRun().getProperty("DANGLE").getStatistics().mean
        dangle0 = ws.getRun().getProperty("DANGLE0").getStatistics().mean
        direct_beam_pix = ws.getRun().getProperty("DIRPIX").getStatistics().mean
        det_distance = ws.getRun().getProperty("SampleDetDis").getStatistics().mean / 1000.0

        peak_pos = peak_position if peak_position is not None else direct_beam_pix
        theta_d = (dangle - dangle0) / 2.0
        theta_d += ((direct_beam_pix - peak_pos) * cls.pixel_width) * 180.0 / math.pi / (2.0 * det_distance)
        return theta_d

    @classmethod
    def scattering_angle_from_data(cls, data_object):
        """
            Compute the scattering angle from a CrossSectionData object, in degrees.
            #TODO: this will go away once we keep the workspace

            @param data_object: CrossSectionData object
        """
        theta_d = (data_object.dangle - data_object.dangle0) / 2.0
        theta_d += ((data_object.dpix - data_object.configuration.peak_position) * cls.pixel_width) * 180.0 / math.pi / (2.0 * data_object.dist_sam_det)
        theta_d += data_object.angle_offset
        return theta_d

    def check_direct_beam(self, ws, peak_position):
        """
            Determine whether this data is a direct beam
        """
        sangle = ws.getRun().getProperty("SANGLE").getStatistics().mean
        theta = self.scattering_angle(ws, peak_position)
        huber_x = ws.getRun().getProperty("HuberX").getStatistics().mean
        return not ((theta > self.tolerance or sangle > self.tolerance) and huber_x < self.huber_x_cut)

    def direct_beam_match(self, scattering, direct_beam, skip_slits=False):
        """
            Verify whether two data sets are compatible.
        """
        if scattering.number == direct_beam.number \
            or ((direct_beam.scattering_angle > self.tolerance \
                 or direct_beam.sangle > self.tolerance) and direct_beam.huber_x < self.huber_x_cut):
            logging.error("Run %s may not be a direct beam", direct_beam.number)

        if math.fabs(scattering.lambda_center-direct_beam.lambda_center) < self.tolerance \
            and (skip_slits or \
            (math.fabs(scattering.slit1_width-direct_beam.slit1_width) < self.tolerance \
            and math.fabs(scattering.slit2_width-direct_beam.slit2_width) < self.tolerance \
            and math.fabs(scattering.slit3_width-direct_beam.slit3_width) < self.tolerance)):
            return True
        return False

    @classmethod
    def get_info(cls, workspace, data_object):
        """
            Retrieve information that is specific to this particular instrument

            @param workspace: Mantid workspace
            @param data_object: CrossSectionData object
        """
        data = workspace.getRun()
        data_object.lambda_center=data['LambdaRequest'].value[0]
        data_object.dangle=data['DANGLE'].value[0]
        data_object.dangle0=data['DANGLE0'].value[0]
        data_object.dpix=data['DIRPIX'].value[0]
        data_object.slit1_width=data['S1HWidth'].value[0]
        data_object.slit2_width=data['S2HWidth'].value[0]
        data_object.slit3_width=data['S3HWidth'].value[0]
        data_object.huber_x=data['HuberX'].getStatistics().mean

        #TODO: these don't exist in the DASLogs
        #data_object.slit1_dist=-data['instrument/aperture1/distance'].value[0]*1000.
        #data_object.slit2_dist=-data['instrument/aperture2/distance'].value[0]*1000.
        #data_object.slit3_dist=-data['instrument/aperture3/distance'].value[0]*1000.

        data_object.sangle=data['SANGLE'].value[0]

        data_object.dist_sam_det=data['SampleDetDis'].value[0]*1e-3
        data_object.dist_mod_det=data['ModeratorSamDis'].value[0]*1e-3+data_object.dist_sam_det
        data_object.dist_mod_mon=data['ModeratorSamDis'].value[0]*1e-3-2.75

        # Get these from instrument
        data_object.pixel_width = float(workspace.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0
        data_object.n_det_size_x = int(workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0]) # 304
        data_object.n_det_size_y = int(workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0]) # 256
        data_object.det_size_x = data_object.n_det_size_x * data_object.pixel_width # horizontal size of detector [m]
        data_object.det_size_y = data_object.n_det_size_y * data_object.pixel_width # vertical size of detector [m]

        # The following active area used to be taken from instrument.DETECTOR_REGION
        data_object.active_area_x = (8, 295)
        data_object.active_area_y = (8, 246)

        # Convert to standard names
        data_object.direct_pixel = data_object.dpix
        data_object.angle_offset = data_object.dangle0

    def process_roi(self, ws):
        """
            Process the ROI information and determine the peak
            range, the low-resolution range, and the background range.
        """
        roi_peak = [0,0]
        roi_low_res = [0,0]
        roi_background = [0,0]

        # Read ROI 1
        roi1_valid = True
        if 'ROI1StartX' in ws.getRun():
            roi1_x0 = ws.getRun()['ROI1StartX'].getStatistics().mean
            roi1_y0 = ws.getRun()['ROI1StartY'].getStatistics().mean
            roi1_x1 = ws.getRun()['ROI1EndX'].getStatistics().mean
            roi1_y1 = ws.getRun()['ROI1EndY'].getStatistics().mean
            if roi1_x1 > roi1_x0:
                peak1 = [int(roi1_x0), int(roi1_x1)]
            else:
                peak1 = [int(roi1_x1), int(roi1_x0)]
            if roi1_y1 > roi1_y0:
                low_res1 = [int(roi1_y0), int(roi1_y1)]
            else:
                low_res1 = [int(roi1_y1), int(roi1_y0)]
            if peak1 == [0,0] and low_res1 == [0,0]:
                roi1_valid = False

            # Read ROI 2
            roi2_valid = True
            roi2_x0 = ws.getRun()['ROI2StartX'].getStatistics().mean
            roi2_y0 = ws.getRun()['ROI2StartY'].getStatistics().mean
            roi2_x1 = ws.getRun()['ROI2EndX'].getStatistics().mean
            roi2_y1 = ws.getRun()['ROI2EndY'].getStatistics().mean
            if roi2_x1 > roi2_x0:
                peak2 = [int(roi2_x0), int(roi2_x1)]
            else:
                peak2 = [int(roi2_x1), int(roi2_x0)]
            if roi2_y1 > roi2_y0:
                low_res2 = [int(roi2_y0), int(roi2_y1)]
            else:
                low_res2 = [int(roi2_y1), int(roi2_y0)]
            if peak2 == [0,0] and low_res2 == [0,0]:
                roi2_valid = False
        else:
            roi1_valid = False
            roi2_valid = False

        # Pick the ROI that describes the reflectivity peak
        if roi1_valid and not roi2_valid:
            roi_peak = peak1
            roi_low_res = low_res1
            roi_background = [0,0]
        elif roi2_valid and not roi1_valid:
            roi_peak = peak2
            roi_low_res = low_res2
            roi_background = [0,0]
        elif roi1_valid and roi2_valid:
            # If ROI 2 is within ROI 1, treat it as the peak,
            # otherwise, use ROI 1
            if peak1[0] >= peak2[0] and peak1[1] <= peak2[1]:
                roi_peak = peak1
                roi_low_res = low_res1
                roi_background = peak2
            elif peak2[0] >= peak1[0] and peak2[1] <= peak1[1]:
                roi_peak = peak2
                roi_low_res = low_res2
                roi_background = peak1
            else:
                roi_peak = peak1
                roi_low_res = low_res1
                roi_background = [0,0]

        return roi_peak, roi_low_res, roi_background

    def integrate_detector(self, ws, specular=True):
        """
            Integrate a workspace along either the main direction (specular=False) or
            the low-resolution direction (specular=True.

            :param ws: Mantid workspace
            :param specular bool: if True, the low-resolution direction is integrated over
        """
        ws_summed = RefRoi(InputWorkspace=ws, IntegrateY=specular,
                           NXPixel=self.n_x_pixel, NYPixel=self.n_y_pixel,
                           ConvertToQ=False,
                           OutputWorkspace="ws_summed")

        integrated = Integration(ws_summed)
        integrated = Transpose(integrated)
        return integrated

    def determine_peak_range(self, integrated_ws, specular=True, max_pixel=230):
        """
            Determine a peak position.
            :param integrated_ws: Mantid workspace, integrated over the low resolution direction.
        """
        x_values = integrated_ws.readX(0)
        y_values = integrated_ws.readY(0)
        e_values = integrated_ws.readE(0)
        try:
            ws_short = CreateWorkspace(DataX=x_values[self.peak_range_offset:max_pixel],
                                       DataY=y_values[self.peak_range_offset:max_pixel],
                                       DataE=e_values[self.peak_range_offset:max_pixel])
            specular_peak, low_res, _ = LRPeakSelection(InputWorkspace=ws_short)
        except:
            logging.error("Peak finding error [specular=%s]: %s" % (specular, sys.exc_value))
            return [0,0]

        if specular:
            peak = [specular_peak[0]+self.peak_range_offset, specular_peak[1]+self.peak_range_offset]
        else:
            # The low-resolution range finder tends to be a bit tight.
            # Broaden it by a third.
            #TODO: Fix the range finder algorithm
            broadening = (low_res[1]-low_res[0])/3.0
            peak = [low_res[0]+self.peak_range_offset-broadening,
                    low_res[1]+self.peak_range_offset+broadening]
        return peak

    @classmethod
    def find_direct_beam(cls):
        """
            Find a direct beam suitable for this data set by looking into
            its data directory.
            :param ws: Mantid workspace
        """
        
        return
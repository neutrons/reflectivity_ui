"""
    Reduction for MR
"""
from __future__ import absolute_import, division, print_function
import sys

# Import mantid according to the application configuration
from . import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *

import numpy as np
import math
import copy
from scipy.optimize import curve_fit
import logging

def _as_ints(a): return [int(a[0]), int(a[1])]

class DataInfo(object):
    n_events_cutoff = 10000

    def __init__(self, ws, cross_section, configuration):
        MRInspectData(Workspace=ws, UseROI=configuration.use_roi,
                      UpdatePeakRange=configuration.update_peak_range,
                      UseROIBck=configuration.use_roi_bck, UseTightBck=configuration.use_tight_bck,
                      BckWidth=int(round(configuration.bck_offset)), HuberXCut=0.0,
                      ForcePeakROI=configuration.force_peak_roi, PeakROI=configuration.peak_roi,
                      ForceLowResPeakROI=configuration.force_low_res_roi, LowResPeakROI=configuration.low_res_roi,
                      ForceBckROI=configuration.force_bck_roi, BckROI=configuration.bck_roi)
                 
        self.cross_section = cross_section
        self.run_number = ws.getRunNumber()

        run_object = ws.getRun()
        self.is_direct_beam = run_object.getProperty("is_direct_beam").value.lower()=='true'
        self.data_type = 0 if self.is_direct_beam else 1
        if ws.getNumberEvents() < self.n_events_cutoff:
            self.data_type = -1

        # Processing options
        # Use the ROI rather than finding the ranges
        self.use_roi = configuration.use_roi
        self.use_roi_actual = run_object.getProperty("use_roi_actual").value.lower()=='true'

        self.calculated_scattering_angle = run_object.getProperty("calculated_scatt_angle").value

        tof_min = run_object.getProperty("tof_range_min").value
        tof_max = run_object.getProperty("tof_range_max").value
        self.tof_range = [tof_min, tof_max]
        
        peak_min = run_object.getProperty("peak_min").value
        peak_max = run_object.getProperty("peak_max").value
        self.peak_range = [peak_min, peak_max]
        self.peak_position = (peak_min+peak_max)/2.0

        background_min = run_object.getProperty("background_min").value
        background_max = run_object.getProperty("background_max").value
        self.background = [background_min, background_max]
        
        low_res_min = run_object.getProperty("low_res_min").value
        low_res_max = run_object.getProperty("low_res_max").value
        self.low_res_range = [low_res_min, low_res_max]

        # Region of interest information
        roi_peak_min = run_object.getProperty("roi_peak_min").value
        roi_peak_max = run_object.getProperty("roi_peak_max").value
        self.roi_peak = [roi_peak_min, roi_peak_max]

        roi_low_res_min = run_object.getProperty("roi_low_res_min").value
        roi_low_res_max = run_object.getProperty("roi_low_res_max").value
        self.roi_low_res = [roi_low_res_min, roi_low_res_max]

        roi_background_min = run_object.getProperty("roi_background_min").value
        roi_background_max = run_object.getProperty("roi_background_max").value
        self.roi_background = [roi_background_min, roi_background_max]


class _DataInfo(object):
    """
        Class to hold the relevant information from a run (scattering or direct beam).
    """
    n_events_cutoff = 10000

    def __init__(self, ws, cross_section, configuration):
        self.cross_section = cross_section
        self.is_direct_beam = False
        self.data_type = 1
        self.peak_position = 0
        self.peak_range = [0,0]
        self.low_res_range = [0,0]
        self.background = [0,0]
        
        # ROI information
        self.roi_peak = [0,0]
        self.roi_low_res = [0,0]
        self.roi_background = [0,0]

        # Options to override the ROI
        self.force_peak_roi = configuration.force_peak_roi
        self.forced_peak_roi = _as_ints(configuration.peak_roi)
        self.force_low_res_roi = configuration.force_low_res_roi
        self.forced_low_res_roi = _as_ints(configuration.low_res_roi)
        self.force_bck_roi = configuration.force_bck_roi
        self.forced_bck_roi = _as_ints(configuration.bck_roi)
        
        # Peak found before fitting for the central position
        self.found_peak = [0,0]
        self.found_low_res = [0,0]

        # Processing options
        self.configuration = configuration
        # Use the ROI rather than finding the ranges
        self.use_roi = configuration.use_roi
        self.use_roi_actual = False

        # Use the 2nd ROI as the background, if available
        self.use_roi_bck = configuration.use_roi_bck

        # Use background as a region on each side of the peak
        self.use_tight_bck = configuration.use_tight_bck
        # Width of the background on each side of the peak
        self.bck_offset = configuration.bck_offset

        # Update the specular peak range after finding the peak
        # within the ROI
        self.update_peak_range = configuration.update_peak_range

        if ws is not None:
            self.run_number = ws.getRunNumber()
            self.tof_range = configuration.instrument.get_tof_range(ws.getRun())
            self.determine_data_type(ws)
        else:
            self.run_number = 0
            self.tof_range = [0, 0]

    def log(self):
        logging.warning("| Run: %s [direct beam: %s]" % (self.run_number, self.is_direct_beam))
        logging.warning("|   Peak position: %s" % self.peak_position)
        logging.warning("|   Reflectivity peak: %s" % str(self.peak_range))
        logging.warning("|   Low-resolution pixel range: %s" % str(self.low_res_range))

    def process_roi(self, ws):
        """
            Process the ROI information and determine the peak
            range, the low-resolution range, and the background range.
        """
        self.roi_peak, self.roi_low_res, self.roi_background = self.configuration.instrument.process_roi(ws)

        # After all this, update the ROI according to reduction options
        if self.force_peak_roi:
            logging.error("Forcing peak ROI: %s", self.forced_peak_roi)
            self.roi_peak = self.forced_peak_roi
        if self.force_low_res_roi:
            logging.error("Forcing low-res ROI: %s", self.forced_low_res_roi)
            self.roi_low_res = self.forced_low_res_roi
        if self.force_bck_roi:
            logging.error("Forcing background ROI: %s", self.forced_bck_roi)
            self.roi_background = self.forced_bck_roi

    @classmethod
    def fit_peak(cls, signal_x, signal_y, peak):
        """
            Find a peak in a data distribution.
            :param signal_x array: array of x values
            :param signal_y array: array of y values
            :param peak array: starting peak position
        """
        def gauss(x, *p):
            A, mu, sigma = p
            return A*np.exp(-(x-mu)**2/(2.*sigma**2))

        p0 = [np.max(signal_y), (peak[1]+peak[0])/2.0, (peak[1]-peak[0])/2.0]
        coeff, _ = curve_fit(gauss, signal_x, signal_y, p0=p0)
        peak_position = coeff[1]
        peak_width = math.fabs(3.0*coeff[2])
        return peak_position, peak_width

    def determine_data_type(self, ws):
        """
            Inspect the data and determine peak locations
            and data type.
        """
        # Skip empty data entries
        if ws.getNumberEvents() < self.n_events_cutoff:
            self.data_type = -1
            logging.info("No data for %s %s" % (self.run_number, self.cross_section))
            return

        # Find reflectivity peak and low resolution ranges
        # Those will be our defaults
        integrated_ws = self.configuration.instrument.integrate_detector(ws, specular=True)
        peak = self.configuration.instrument.determine_peak_range(integrated_ws, specular=True)
        self.found_peak = copy.copy(peak)
        logging.info("Run %s [%s]: Peak found %s" % (self.run_number, self.cross_section, peak))
        signal_y = integrated_ws.readY(0)
        signal_x = range(len(signal_y))

        integrated_ws = self.configuration.instrument.integrate_detector(ws, specular=False)
        low_res = self.configuration.instrument.determine_peak_range(integrated_ws, specular=False)
        logging.info("Run %s [%s]: Low-res found %s" % (self.run_number, self.cross_section, low_res))
        self.found_low_res = low_res
        bck_range = None
        
        # Process the ROI information
        self.process_roi(ws)

        # Keep track of whether we actually used the ROI
        self.use_roi_actual = False
        
        if self.use_roi and not self.roi_peak == [0,0]:
            peak = copy.copy(self.roi_peak)
            if not self.roi_low_res == [0,0]:
                low_res = copy.copy(self.roi_low_res)
            if not self.roi_background == [0,0]:
                bck_range = copy.copy(self.roi_background)
            logging.info("Using ROI peak range: [%s %s]", peak[0], peak[1])
            self.use_roi_actual = True

        # Determine reflectivity peak position (center)
        signal_y_crop = signal_y[peak[0]:peak[1]+1]
        signal_x_crop = signal_x[peak[0]:peak[1]+1]

        peak_position = (peak[1]+peak[0])/2.0
        peak_width = (peak[1]-peak[0])/2.0
        try:
            # Try to find the peak position within the peak range we found
            peak_position, peak_width = self.fit_peak(signal_x_crop, signal_y_crop, peak)
            # If we are more than two sigmas away from the middle of the range,
            # there's clearly a problem.
            if np.abs(peak_position - (peak[1]+peak[0])/2.0)  > np.abs(peak[1]-peak[0]):
                logging.error("Found peak position outside of given range [x=%s], switching to full detector" % peak_position)
                peak_position = (peak[1]+peak[0])/2.0
                peak_width = (peak[1]-peak[0])/2.0
                raise RuntimeError("Bad peak position")
        except:
            # If we can't find a peak, try fitting over the full detector.
            # If we do find a peak, then update the ranges rather than using
            # what we currently have (which is probably given by the ROI).
            logging.warning("Run %s [%s]: Could not fit a peak in the supplied peak range" % (self.run_number, self.cross_section))
            logging.warning(sys.exc_value)
            try:
                # Define a good default that is wide enough for the fit to work
                default_width = (self.found_peak[1]-self.found_peak[0])/2.0
                default_width = max(default_width, 5.0)
                default_center = (self.found_peak[1]+self.found_peak[0])/2.0
                default_peak = [default_center-default_width, default_center+default_width]
                peak_position, peak_width = self.fit_peak(signal_x, signal_y, default_peak)
                peak = [math.floor(peak_position-peak_width), math.floor(peak_position+peak_width)]
                #low_res = [5, self.n_x_pixel-5]
                low_res = self.found_low_res
                self.use_roi_actual = False
                logging.warning("Run %s [%s]: Peak not in supplied range! Found peak: %s low: %s" % (self.run_number, self.cross_section, peak, low_res))
                logging.warning("Run %s [%s]: Peak position: %s  Peak width: %s" % (self.run_number, self.cross_section, peak_position, peak_width))
            except:
                logging.warning(sys.exc_value)
                logging.error("Run %s [%s]: Could not use Gaussian fit to determine peak position over whole detector" % (self.run_number, self.cross_section))

        # Update the specular peak range if needed
        if self.update_peak_range:
            peak[0] = math.floor(peak_position-peak_width)
            peak[1] = math.ceil(peak_position+peak_width)
            logging.info("Updating peak range to: [%s %s]", peak[0], peak[1])
            self.use_roi_actual = False

        # Store the information we found
        self.peak_position = peak_position
        self.peak_range = [int(peak[0]), int(peak[1])]
        self.low_res_range = [int(low_res[0]), int(low_res[1])]

        if not self.use_roi_bck or bck_range is None:
            if self.use_tight_bck:
                self.background = [self.peak_range[0]-self.bck_offset, self.peak_range[1]+self.bck_offset]
            else:
                self.background = [4, self.peak_range[0]-30]
        else:
            self.background = [int(bck_range[0]), int(bck_range[1])]

        # Computed scattering angle
        self.scattering_angle = self.configuration.instrument.scattering_angle(ws, peak_position)
        
        # Determine whether we have a direct beam
        self.is_direct_beam = self.configuration.instrument.check_direct_beam(ws, peak_position)

        # Convenient data type
        self.data_type = 0 if self.is_direct_beam else 1

        # Write to logs
        self.log()


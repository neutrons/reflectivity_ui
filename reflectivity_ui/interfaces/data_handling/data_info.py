"""
    Meta-data information for MR reduction
"""
#pylint: disable=too-few-public-methods, wrong-import-position, too-many-instance-attributes, wrong-import-order
from __future__ import absolute_import, division, print_function
import sys
import time
import logging
import math
import copy

import numpy as np
import scipy.optimize as opt

# Import mantid according to the application configuration
from . import ApplicationConfiguration
sys.path.insert(0, ApplicationConfiguration().mantid_path)
import mantid.simpleapi as api


NX_PIXELS = 304
NY_PIXELS = 256


class DataInfo(object):
    """
        Class to hold the relevant information from a run (scattering or direct beam).
    """
    peak_range_offset = 0
    tolerance = 0.02

    def __init__(self, ws, cross_section, configuration):
        self.cross_section = cross_section
        self.run_number = ws.getRunNumber()
        self.is_direct_beam = False
        self.data_type = 1
        self.peak_position = 0
        self.peak_range = [0,0]
        self.low_res_range = [0,0]
        self.background = [0,0]
        self.n_events_cutoff = 10000

        # ROI information
        self.roi_peak = [0,0]
        self.roi_low_res = [0,0]
        self.roi_background = [0,0]

        # Options to override the ROI
        self.force_peak_roi = configuration.force_peak_roi
        self.forced_peak_roi = configuration.peak_roi
        self.force_low_res_roi = configuration.force_low_res_roi
        self.forced_low_res_roi = configuration.low_res_roi
        self.force_bck_roi = configuration.force_bck_roi
        self.forced_bck_roi = configuration.bck_roi

        # Peak found before fitting for the central position
        self.found_peak = [0,0]
        self.found_low_res = [0,0]

        # Processing options
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

        self.tof_range = self.get_tof_range(ws)
        self.calculated_scattering_angle = 0.0
        self.theta_d = 0.0
        t_0 = time.time()
        self.determine_data_type(ws)
        logging.info("INSPECT: %s sec" % (time.time()-t_0))

    def get_tof_range(self, ws):
        """
            Determine TOF range from the data
            :param workspace ws: workspace to work with
        """
        run_object = ws.getRun()
        sample_detector_distance = run_object['SampleDetDis'].getStatistics().mean
        source_sample_distance = run_object['ModeratorSamDis'].getStatistics().mean
        # Check units
        if not run_object['SampleDetDis'].units in ['m', 'meter']:
            sample_detector_distance /= 1000.0
        if not run_object['ModeratorSamDis'].units in ['m', 'meter']:
            source_sample_distance /= 1000.0

        source_detector_distance = source_sample_distance + sample_detector_distance

        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        wl = run_object.getProperty('LambdaRequest').value[0]
        chopper_speed = run_object.getProperty('SpeedRequest1').value[0]
        wl_offset = 0
        cst = source_detector_distance / h * m
        tof_min = cst * (wl + wl_offset * 60.0 / chopper_speed - 1.4 * 60.0 / chopper_speed) * 1e-4
        tof_max = cst * (wl + wl_offset * 60.0 / chopper_speed + 1.4 * 60.0 / chopper_speed) * 1e-4

        self.tof_range = [tof_min, tof_max]
        return [tof_min, tof_max]

    def process_roi(self, ws):
        """
            Process the ROI information and determine the peak
            range, the low-resolution range, and the background range.

            Starting in June 2018, with the DAS upgrade, the ROIs are
            specified with a start/width rather than start/stop.

            :param workspace ws: workspace to work with
        """
        roi_peak = [0,0]
        roi_low_res = [0,0]
        roi_background = [0,0]

        # Read ROI 1
        roi1_valid = True
        if 'ROI1StartX' in ws.getRun():
            roi1_x0 = ws.getRun()['ROI1StartX'].getStatistics().mean
            roi1_y0 = ws.getRun()['ROI1StartY'].getStatistics().mean
            if 'ROI1SizeX' in ws.getRun():
                size_x = ws.getRun()['ROI1SizeX'].getStatistics().mean
                size_y = ws.getRun()['ROI1SizeY'].getStatistics().mean
                roi1_x1 = roi1_x0 + size_x
                roi1_y1 = roi1_y0 + size_y
            else:
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
            if 'ROI2StartX' in ws.getRun():
                roi2_valid = True
                roi2_x0 = ws.getRun()['ROI2StartX'].getStatistics().mean
                roi2_y0 = ws.getRun()['ROI2StartY'].getStatistics().mean
                if 'ROI2SizeX' in ws.getRun():
                    size_x = ws.getRun()['ROI2SizeX'].getStatistics().mean
                    size_y = ws.getRun()['ROI2SizeY'].getStatistics().mean
                    roi2_x1 = roi2_x0 + size_x
                    roi2_y1 = roi2_y0 + size_y
                else:
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

        # After all this, update the ROI according to reduction options
        self.roi_peak = roi_peak
        self.roi_low_res = roi_low_res
        self.roi_background = roi_background

    def determine_data_type(self, ws):
        """
            Inspect the data and determine peak locations
            and data type.
            :param workspace ws: Workspace to inspect
        """
        # Skip empty data entries
        if ws.getNumberEvents() < self.n_events_cutoff:
            self.data_type = -1
            logging.info("No data for %s %s" % (self.run_number, self.cross_section))
            return

        # Find reflectivity peak and low resolution ranges
        fitter = Fitter(ws, True)
        peak, low_res = fitter.fit_2d_peak()

        if self.use_tight_bck:
            bck_range = [int(max(0.0, peak[0]-self.bck_offset)), int(min(NX_PIXELS, peak[1]+self.bck_offset))]
        else:
            bck_range = [int(max(0.0, peak[0]-2*self.bck_offset)), int(max(0.0, peak[0]-self.bck_offset))]
        self.found_peak = copy.copy(peak)
        self.found_low_res = copy.copy(low_res)
        logging.info("Run %s [%s]: Peak found %s" % (self.run_number, self.cross_section, peak))
        logging.info("Run %s [%s]: Low-res found %s" %(self.run_number, self.cross_section, str(low_res)))

        # Process the ROI information
        try:
            self.process_roi(ws)
        except:
            logging.info("Could not process ROI\n%s" % sys.exc_info()[1])

        # Keep track of whether we actually used the ROI
        self.use_roi_actual = False

        # If we were asked to use the ROI but no peak is in it, use the peak we found
        # If we were asked to use the ROI and there's a peak in it, use the ROI
        if self.use_roi and not self.update_peak_range and not self.roi_peak == [0,0]:
            logging.info("Using ROI peak range: [%s %s]" % (self.roi_peak[0], self.roi_peak[1]))
            self.use_roi_actual = True
            peak = copy.copy(self.roi_peak)
            if not self.roi_low_res == [0,0]:
                low_res = copy.copy(self.roi_low_res)
            if not self.roi_background == [0,0]:
                bck_range = copy.copy(self.roi_background)
        elif self.use_roi and self.update_peak_range and not self.roi_peak == [0,0]:
            logging.info("Using fit peak range: [%s %s]" % (peak[0], peak[1]))
            if not self.roi_background == [0,0]:
                bck_range = copy.copy(self.roi_background)

        # Store the information we found
        self.peak_position = (peak[1]+peak[0])/2.0
        self.peak_range = [int(max(0, peak[0])), int(min(peak[1], NX_PIXELS))]
        self.low_res_range = [int(max(0, low_res[0])), int(min(low_res[1], NY_PIXELS))]
        self.background = [int(max(0, bck_range[0])), int(min(bck_range[1], NY_PIXELS))]

        # Computed scattering angle
        self.calculated_scattering_angle = api.MRGetTheta(ws, SpecularPixel=self.peak_position)
        self.calculated_scattering_angle *= 180.0 / math.pi

        # Determine whether we have a direct beam
        run_object = ws.getRun()
        try:
            self.is_direct_beam = run_object.getProperty("data_type").value[0] == 1
            self.data_type = 0 if self.is_direct_beam else 1
        except:
            self.is_direct_beam = False
            self.data_type = 1


class DataInfo_(object):
    """
        Class to provide a convenient interface to the meta-data extracted
        by MRInspectData.
    """
    # Cutoff below which we can't call a data set a direct beam
    #n_events_cutoff = 2000

    def __init__(self, ws, cross_section, configuration):
        t_0 = time.time()
        api.MRInspectData(Workspace=ws, UseROI=configuration.use_roi,
                          UpdatePeakRange=configuration.update_peak_range,
                          UseROIBck=configuration.use_roi_bck, UseTightBck=configuration.use_tight_bck,
                          BckWidth=int(round(configuration.bck_offset)),
                          ForcePeakROI=configuration.force_peak_roi, PeakROI=configuration.peak_roi,
                          ForceLowResPeakROI=configuration.force_low_res_roi, LowResPeakROI=configuration.low_res_roi,
                          ForceBckROI=configuration.force_bck_roi, BckROI=configuration.bck_roi)
        t_1 = time.time()
        self.cross_section = cross_section
        self.run_number = ws.getRunNumber()

        run_object = ws.getRun()
        try:
            self.is_direct_beam = run_object.getProperty("data_type").value[0] == 1
            self.data_type = 0 if self.is_direct_beam else 1
        except:
            self.is_direct_beam = False
            self.data_type = 1

        # Processing options
        # Use the ROI rather than finding the ranges
        self.use_roi = configuration.use_roi
        self.use_roi_actual = run_object.getProperty("use_roi_actual").value.lower() == 'true'

        self.calculated_scattering_angle = run_object.getProperty("calculated_scatt_angle").value

        tof_min = run_object.getProperty("tof_range_min").value
        tof_max = run_object.getProperty("tof_range_max").value
        self.tof_range = [tof_min, tof_max]

        # Region of interest information
        roi_peak_min = run_object.getProperty("roi_peak_min").value
        roi_peak_max = run_object.getProperty("roi_peak_max").value
        self.roi_peak = [roi_peak_min, roi_peak_max]

        improved_peaks = True
        if improved_peaks:
            fitter = Fitter(ws, False)
            peak_min = run_object.getProperty("peak_min").value
            peak_max = run_object.getProperty("peak_max").value
            #[peak_min, peak_max], [low_res_min, low_res_max] = fitter.fit_2d_peak()
            [low_res_min, low_res_max] = fitter.fit_beam_width()
            if np.abs(peak_max-peak_min)<=1:
                    peak_min = peak_min-2
                    peak_max = peak_max+2
        else:
            peak_min = run_object.getProperty("peak_min").value
            peak_max = run_object.getProperty("peak_max").value
            low_res_min = run_object.getProperty("low_res_min").value
            low_res_max = run_object.getProperty("low_res_max").value

        self.peak_range = [peak_min, peak_max]
        self.peak_position = (peak_min+peak_max)/2.0
        self.low_res_range = [low_res_min, low_res_max]

        bck_max = min(20, self.peak_position)
        self.background = [max(0, bck_max-10), bck_max]

        roi_low_res_min = run_object.getProperty("roi_low_res_min").value
        roi_low_res_max = run_object.getProperty("roi_low_res_max").value
        self.roi_low_res = [roi_low_res_min, roi_low_res_max]

        roi_background_min = run_object.getProperty("roi_background_min").value
        roi_background_max = run_object.getProperty("roi_background_max").value
        self.roi_background = [roi_background_min, roi_background_max]
        logging.info("DataInfo: %s sec [MRInspectData: %s sec]", (time.time() - t_0), (t_1 - t_0))


def chi2(data, model):
    """ Returns the chi^2 for a data set and model pair """
    err = np.fabs(data)
    err[err<=0] = 1
    return np.sum((data - model)**2 / err) / len(data)


class Fitter(object):
    """
        Peak finder for MR data
    """
    DEAD_PIXELS = 10
    DEFAULT_PEAK_WIDTH = 3

    def __init__(self, workspace, prepare_plot_data=False):
        self.workspace = workspace
        self.prepare_plot_data = prepare_plot_data
        self._prepare_data()
        logging.info("Numpy version: %s" % np.__version__)

    def _prepare_data(self):
        """
            Read in the data and create arrays for fitting
        """
        # Prepare data to fit
        self.n_x = int(self.workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0])
        self.n_y = int(self.workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0])
        self.dirpix = self.workspace.getRun()['DIRPIX'].value[0]

        _integrated = api.Integration(InputWorkspace=self.workspace)
        signal = _integrated.extractY()
        self.z=np.reshape(signal, (self.n_x, self.n_y))
        self.x = np.arange(0, self.n_x)[self.DEAD_PIXELS:-self.DEAD_PIXELS]
        self.y = np.arange(0, self.n_y)[self.DEAD_PIXELS:-self.DEAD_PIXELS]

        # 1D data x/y vs counts
        self.x_vs_counts = np.sum(self.z, 1)[self.DEAD_PIXELS:-self.DEAD_PIXELS]
        self.y_vs_counts = np.sum(self.z, 0)[self.DEAD_PIXELS:-self.DEAD_PIXELS]
        self.x_err = np.sqrt(self.x_vs_counts)
        self.x_err[self.x_err<1] = 1

        # Use the highest data point as a starting point for a simple Gaussian fit
        self.center_x = np.argmax(self.x_vs_counts)
        self.center_y = np.argmax(self.y_vs_counts)

        self.guess_x = self.center_x
        self.guess_wx = 6
        self.guess_y = self.center_y
        self.guess_wy = 50
        self.guess_chi2 = np.inf

        # Plots [optional]
        self.plot_list = []
        self.plot_labels = []
        if self.prepare_plot_data:
            self.plot_list = [[self.x, self.x_vs_counts],]
            self.plot_labels = ['Data',]

    def _perform_beam_fit(self, y_d, derivative, derivative_err, y_r=None, signal_r=None, gaussian_first=False):
        if gaussian_first:
            _running_err = np.sqrt(signal_r)
            _gauss, _ = opt.curve_fit(self.gaussian_1d, y_r,
                                      signal_r, p0=[np.max(signal_r), 140, 50, 0], sigma=_running_err)
            p0 = [np.max(derivative), _gauss[1], 2.0*_gauss[2], 5, 0]
        else:
            p0 = [np.max(derivative), 140, 60, 5, 0]

        #p = A, center_x, width_x, edge_width, background
        _coef, _ = opt.curve_fit(self.peak_derivative, y_d, derivative, p0=p0, sigma=derivative_err)
        return _coef

    def fit_beam_width(self):
        """
            Fit the data distribution in y and get its range.
        """
        _integral = [np.sum(self.y_vs_counts[:i]) for i in range(len(self.y_vs_counts))]
        _running = 0.1*np.convolve(self.y_vs_counts, np.ones(10), mode='valid')
        _deriv = np.asarray([_running[i+1]-_running[i] for i in range(len(_running)-1)])
        _deriv_err = np.sqrt(_running)[:-1]
        _deriv_err[_deriv_err<1] = 1
        _y = self.y[5:-5]

        _coef = self._perform_beam_fit(_y, _deriv, _deriv_err, gaussian_first=False)
        peak_min = _coef[1] - np.abs(_coef[2])/2.0 - 2.0 * np.abs(_coef[3])
        peak_max = _coef[1] + np.abs(_coef[2])/2.0 + 2.0 * np.abs(_coef[3])
        if peak_max - peak_min < 10:
            logging.error("Low statisting: trying again")
            _y_running = self.y[5:-4]
            _coef = self._perform_beam_fit(_y, _deriv, _deriv_err, _y_running, _running, gaussian_first=True)

        self.guess_y = _coef[1]
        self.guess_wy = (peak_max - peak_min) / 2.0
        return [peak_min, peak_max]

    def get_roi(self, region):
        """
            Select are region of interest and prepare the data for fitting.
            :param region: Length 2 list of min/max pixels defining the ROI
        """
        _roi = [x_i>region[0] and x_i<=region[1] for x_i in self.x]
        x_roi = self.x[_roi]
        counts_roi = self.x_vs_counts[_roi]
        err_roi = np.sqrt(np.fabs(counts_roi))
        err_roi[err_roi<1] = 1

        return x_roi, counts_roi, err_roi

    def _scan_peaks(self):
        """
            Perform a quick scan of the count distribution in x to find obvious peaks.
            We first convolute the distribution with a step function to get rid of
            noise, then we compute the first derivative. We identify a peak when the
            difference in the slop between two consecutive points is greater than
            half the total counts.
        """
        found_peaks = []
        # Derivative
        _convo_narrow = np.zeros(len(self.x_vs_counts))
        _width = 2
        total_counts = np.sum(self.x_vs_counts)

        for i, _ in enumerate(_convo_narrow):
            _step = [1 if abs(i-j)<=_width else 0 for j in range(len(self.x_vs_counts))]
            _convo_narrow[i] = np.sum(_step*self.x_vs_counts)
        _deriv = [ _convo_narrow[i+1] - _convo_narrow[i] for i in range(len(_convo_narrow)-1) ]

        _up = None
        _i_value = 0
        for i in range(len(_deriv)-1):
            if i==0:
                _up = _deriv[i+1] > _deriv[i]
                _i_value = i + 1
            elif _up:
                if _deriv[i+1] > _deriv[i]:
                    _i_value = i + 1
                else:
                    _up = False
            else:
                if _deriv[i+1] > _deriv[i]:
                    if _deriv[_i_value] - _deriv[i+1] > np.sqrt(total_counts)/2 \
                    and _i_value>self.DEAD_PIXELS and _i_value<self.n_x-self.DEAD_PIXELS:
                        found_peaks.append(int((_i_value + i) / 2.0))
                    _up = True
                    _i_value = i + 1

        # Too many peaks is not good, because it means that we haven't correctly
        # identified the reflected beam and the direct beam.
        if len(found_peaks)>2:
            return []

        if found_peaks:
            self.guess_x = self.x[max(1, found_peaks[0] - self.DEFAULT_PEAK_WIDTH)]
            self.guess_wx = self.x[min(self.n_x-2, found_peaks[0] + self.DEFAULT_PEAK_WIDTH)]

        peaks = [self.x[i] for i in found_peaks]
        return peaks

    def _fit_gaussian(self):
        """
            Fit a simple Gaussian and constant background
        """
        if self.peaks:
            center_x = self.peaks[0]
        else:
            center_x = self.center_x

        # Scale, mu_x, sigma_x, mu_y, sigma_y, background
        p0 = [np.max(self.x_vs_counts), center_x, 5, 0]
        try:
            gauss_coef, _ = opt.curve_fit(self.gaussian_1d,
                                          self.x,
                                          self.x_vs_counts, p0=p0,
                                          sigma=self.x_err)
        except:
            logging.info("Could not fit simple Gaussian")
            gauss_coef = p0

        # Keep track of the result
        theory = self.gaussian_1d(self.x, *gauss_coef)
        _chi2 = chi2(theory, self.x_vs_counts)

        # Fitting a Gaussian tends to give a narrower peak than we
        # really need, so we're multiplying the width by two.
        if _chi2 < self.guess_chi2:
            self.guess_x = gauss_coef[1]
            self.guess_wx = 2.0 * gauss_coef[2]
            self.guess_chi2 = _chi2

        if self.prepare_plot_data:
            self.plot_list.append([self.x, theory])
            self.plot_labels.append('Gaussian')
            logging.info("Chi2[Gaussian] = %s" % _chi2)
            logging.info("    %g +- %g" % (gauss_coef[1], gauss_coef[2]))

    def _fit_gaussian_and_poly(self):
        """
            Fit a Gaussian with a polynomial background. First fit the background,
            then keep it constant and add a Gaussian.
        """
        if self.peaks:
            center_x = self.peaks[0]
        else:
            center_x = self.center_x

        try:
            poly_bck_coef, _ = opt.curve_fit(self.poly_bck, self.x, self.x_vs_counts,
                                             p0=[np.max(self.x_vs_counts), 0, 0, center_x, 0], sigma=self.x_err)
        except:
            logging.info("Could not fit polynomial background")
            poly_bck_coef = [0, 0, 0, self.center_x, 0]
        theory = self.poly_bck(self.x, *poly_bck_coef)

        if self.prepare_plot_data:
            _chi2 = chi2(theory, self.x_vs_counts)
            self.plot_list.append([self.x, theory])
            self.plot_labels.append('Polynomial')
            logging.info("Chi2[Polynomial] = %g" % _chi2)

        # Now fit a Gaussian + background
        # A, mu_x, sigma_x, background
        self.poly_bck_coef = poly_bck_coef
        coef = [np.max(self.x_vs_counts), self.center_x, 5, 0]
        try:
            coef, _ = opt.curve_fit(self.gaussian_and_fixed_poly_bck, self.x, self.x_vs_counts,
                                    p0=coef, sigma=self.x_err)
        except:
            logging.info("Could not fit Gaussian + polynomial")
        theory = self.gaussian_and_fixed_poly_bck(self.x, *coef)
        _chi2 = chi2(theory, self.x_vs_counts)

        # Fitting a Gaussian tends to give a narrower peak than we
        # really need, so we're multiplying the width by two.
        if _chi2 < self.guess_chi2:
            self.guess_x = coef[1]
            self.guess_wx = 2.0 * coef[2]
            self.guess_chi2 = _chi2

        if self.prepare_plot_data:
            self.plot_list.append([self.x, theory])
            self.plot_labels.append('Gaussian + polynomial')
            logging.info("Chi2[Gaussian + polynomial] = %g" % _chi2)
            logging.info("    %g +- %g" % (coef[1], coef[2]))

    def _fit_lorentz_2d(self, peak=True):
        """
            Fit a Lorentzian peak, usually to fit the direct beam
        """
        # Find good starting point according to whether we are fitting a
        # reflected peak or a direct beam peak.
        dirpix= self.center_x
        if self.peaks and len(self.peaks) < 3:
            dirpix = self.peaks[0]
        if not peak:
            if len(self.peaks) > 1 and len(self.peaks) < 3:
                dirpix = self.peaks[1]
            else:
                dirpix = self.dirpix

        # Scale, mu_x, fwhm, mu_y, sigma_y, background
        p0 = [np.max(self.x_vs_counts), dirpix, 10, 0]
        try:
            lorentz_coef, _ = opt.curve_fit(self.lorentzian,
                                            self.x,
                                            self.x_vs_counts, p0=p0,
                                            sigma=self.x_err)
        except:
            logging.info("Could not fit Lorentzian")
            lorentz_coef = p0

        # Keep track of the result
        theory = self.lorentzian(self.x, *lorentz_coef)
        _chi2 = chi2(theory, self.x_vs_counts)

        if self.prepare_plot_data:
            self.plot_list.append([self.x, theory])
            self.plot_labels.append('Lorentz 2D')
            logging.info("Chi2[Lorentz 2D] = %s" % _chi2)
            logging.info("    %g +- %g" % (lorentz_coef[1], lorentz_coef[2]))
        return lorentz_coef

    def _gaussian_and_lorentzian(self, region):
        """
            Fit a Gaussian on top of a lorentzian background.
            First fit the Lorentzian, then keep it constant to fit a Gaussian on top.
            :param list region: Length 2 list defining the x region to fit in
        """
        if len(self.peaks) > 1:
            center_x = self.peaks[0]
        elif len(self.peaks) == 1:
            center_x = self.peaks[0]
        else:
            center_x = self.center_x

        # Fit the Lorentzian background first
        self.lorentz_coef = self._fit_lorentz_2d(peak=False)

        # Extract the region we want to fit over
        code_roi, data_to_fit_roi, err_roi = self.get_roi(region)

        #A, mu_x, sigma_x, poly_a, poly_b, poly_c, center, background
        p0 = [np.max(data_to_fit_roi), center_x, 5, 0, 0, 0, center_x, 0]
        try:
            lorentz_coef, _ = opt.curve_fit(self.gaussian_and_fixed_lorentzian,
                                            code_roi,
                                            data_to_fit_roi, p0=p0,
                                            sigma=err_roi)
        except:
            logging.info("Could not fit G+L")
            lorentz_coef = p0

        logging.info("G+L params: %s" % str(lorentz_coef))
        # Keep track of the result
        theory = self.gaussian_and_fixed_lorentzian(self.x, *lorentz_coef)
        _chi2 = chi2(theory, self.x_vs_counts)

        # If we decided to fit two peaks, we should take the results regardless
        # of goodness of fit because the models are imprecise.
        # Nonetheless, log an entry if the chi^2 is larger
        if _chi2 > self.guess_chi2:
            logging.info("Fitting with two peaks resulted in a larger chi^2: %g > %g" % (_chi2, self.guess_chi2))

        # Unless we have a crazy peak
        if lorentz_coef[1] > self.peaks[0]-10 and lorentz_coef[1] < self.peaks[0]+10:
            # Fitting a Gaussian tends to give a narrower peak than we
            # really need, so we're multiplying the width by two.
            self.guess_x = lorentz_coef[1]
            self.guess_wx = 2.0 * lorentz_coef[2]
            self.guess_chi2 = _chi2

        if self.prepare_plot_data:
            self.plot_list.append([self.x, theory])
            self.plot_labels.append('G + Lorentz 2D')
            logging.info("Chi2[G + Lorentz] = %s" % _chi2)
            logging.info("    %g +- %g" % (lorentz_coef[1], lorentz_coef[2]))
        return lorentz_coef

    def fit_2d_peak(self, region=None):
        """
            Fit a 2D Gaussian peak
            :param region: region of interest for the reflected peak
        """
        self.peaks = self._scan_peaks()
        logging.info("Peaks (rough scan): %s" % self.peaks)

        # Gaussian fit
        self._fit_gaussian()

        # Fit a polynomial background, as a starting point to fitting signal + background
        self._fit_gaussian_and_poly()

        if len(self.peaks) > 1:
            if region is None:
                region = [self.peaks[0]-20, self.peaks[0]+20]
            self._gaussian_and_lorentzian(region)

        # Fit the beam width (low-res direction)
        self.fit_beam_width()

        # Package the best results
        x_min = max(0, int(self.guess_x-np.fabs(self.guess_wx)))
        x_max = min(self.n_x-1, int(self.guess_x+np.fabs(self.guess_wx)))
        y_min = max(0, int(self.guess_y-np.fabs(self.guess_wy)))
        y_max = min(self.n_y-1, int(self.guess_y+np.fabs(self.guess_wy)))

        return [x_min, x_max], [y_min, y_max]

    # Fit function definitions #####################################################
    def poly_bck(self, value, *p):
        """
            Polynomial function for background fit

            f = a + b*(x-center) + c*(x-center)**2 + bck

            where bck is a minimum threshold that is zero when the polynomial
            has a value greater than it.
        """
        poly_a, poly_b, poly_c, center, background = p
        values = poly_a + poly_b*(value-center) + poly_c*(value-center)**2 + background
        return values

    def gaussian_and_poly_bck(self, value, *p):
        """
            Function for a polynomial + Gaussian signal
        """
        A, mu_x, sigma_x, poly_a, poly_b, poly_c, center, background = p
        poly_coef = [poly_a, poly_b, poly_c, center, background]
        values = self.poly_bck(value, *poly_coef)
        gauss_coef = [A, mu_x, sigma_x, 0]
        values += self.gaussian_1d(value, *gauss_coef)
        return values

    def gaussian_and_fixed_poly_bck(self, value, *p):
        """
            Use result of bck fit and add a Gaussian
        """
        values = self.poly_bck(value, *self.poly_bck_coef)
        values += self.gaussian_1d(value, *p)
        return values

    def lorentzian(self, value, *p):
        """
            Peak function in 2D. The main axis (x) is a Lorentzian and the other axis (y) is a Gaussian.
        """
        A, mu_x, sigma_x, background = p
        values =  abs(A)/(1+((value-mu_x)/sigma_x)**2) + abs(background)
        return values

    def gaussian_and_fixed_lorentzian(self, value, *p):
        """
            Gaussian and polynomial on top of a fixed Lorentzian.
        """
        values = self.lorentzian(value, *self.lorentz_coef)
        values += self.gaussian_and_poly_bck(value, *p)
        return values

    def gaussian_1d(self, value, *p):
        """
            1D Gaussian
        """
        A, center_x, width_x, background = p
        A = np.abs(A)
        values = A*np.exp(-(value-center_x)**2/(2.*width_x**2))
        values += background
        return values

    def peak_derivative(self, value, *p):
        """
            Double Gaussian to fit the first derivative of a plateau/peak.
        """
        A, center_x, width_x, edge_width, background = p
        mu_right = center_x + width_x / 2.0
        mu_left = center_x - width_x / 2.0
        A = np.abs(A)
        values = A*np.exp(-(value-mu_left)**2/(2.*edge_width**2)) - A*np.exp(-(value-mu_right)**2/(2.*edge_width**2))
        values += background
        return values

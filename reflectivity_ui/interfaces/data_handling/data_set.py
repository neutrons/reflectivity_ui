"""
    Loader for event nexus files.
    Uses Mantid Framework
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long, multiple-statements
from __future__ import absolute_import, division, print_function
import sys
import os
import logging
from collections import OrderedDict
import h5py
import numpy as np
import copy

# Import mantid according to the application configuration
from . import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *
# Set Mantid logging level to warnings
ConfigService.setConsoleLogLevel(4)

from .data_info import DataInfo

### Parameters needed for some calculations.

ANALYZER_IN=(0., 100.) # position and maximum deviation of analyzer in it's working position
POLARIZER_IN=(-348., 50.) # position and maximum deviation of polarizer in it's working position
SUPERMIRROR_IN=(19.125, 10.) # position and maximum deviation of the supermirror translation
POLY_CORR_PARAMS=[-4.74152261e-05,-4.62469580e-05, 1.25995446e-02, 2.13654008e-02,
                  1.02334517e+01] # parameters used in polynomial detector sensitivity correction

XSECT_MAPPING = {u'entry-Off_Off': u'++',
                 u'entry-On_On': u'--',
                 u'entry-Off_On': u'+-',
                 u'entry-On_Off': u'-+',
                }

H_OVER_M_NEUTRON = 3.956034e-7 # h/m_n [m^2/s]

def getIxyt(nxs_data):
    """
        Return [x, y, TOF] array
        @param nxs_data: Mantid workspace
    """
    _tof_axis = nxs_data.readX(0)[:].copy()
    nbr_tof = len(_tof_axis)

    sz_y_axis = int(nxs_data.getInstrument().getNumberParameter("number-of-y-pixels")[0]) #256
    sz_x_axis = int(nxs_data.getInstrument().getNumberParameter("number-of-x-pixels")[0]) #304

    _y_axis = np.zeros((sz_x_axis, sz_y_axis, nbr_tof-1))
    #_y_error_axis = np.zeros((sz_x_axis, sz_y_axis, nbr_tof-1))

    for x in range(sz_x_axis):
        for y in range(sz_y_axis):
            _index = int(sz_y_axis*x+y)
            _tmp_data = nxs_data.readY(_index)[:]
            _y_axis[x,y,:] = _tmp_data

    return _y_axis

class NexusData(object):
    """
        Read a nexus file with multiple cross-section data.
    """
    def __init__(self, file_path, configuration):
        self.file_path = file_path
        self.number = 0
        self.configuration = configuration
        self.cross_sections = {}

    @property
    def nbytes(self):
        total_size = 0
        for d in self.cross_sections.keys():
            total_size += self.cross_sections[d].nbytes
        return total_size

    def get_q_range(self):
        """
            Return the Q range for the cross-sections
        """
        q_min = None
        q_max = None
        for xs in self.cross_sections:
            if self.cross_sections[xs].q is not None:
                if q_min is None:
                    q_min = self.cross_sections[xs].q.min()
                    q_max = self.cross_sections[xs].q.max()
                else:
                    q_min = min(q_min, self.cross_sections[xs].q.min())
                    q_max = max(q_max, self.cross_sections[xs].q.max())
        return q_min, q_max

    def set_parameter(self, param, value):
        """
            Loop through the cross-section data sets and update
            a parameter.
        """
        has_changed = False
        try:
            for xs in self.cross_sections:
                if hasattr(self.cross_sections[xs].configuration, param):
                    if not getattr(self.cross_sections[xs].configuration, param) == value:
                        setattr(self.cross_sections[xs].configuration, param, value)
                        has_changed = True
        except:
            logging.error("Could not set parameter %s %s\n  %s", param, value, sys.exc_value)
        return has_changed

    def calculate_reflectivity(self, direct_beam=None, configuration=None):
        """
            Loop through the cross-section data sets and update
            the reflectivity.
        """
        has_errors = False
        detailed_msg = ""
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].reflectivity(direct_beam=direct_beam, configuration=configuration)
            except:
                has_errors = True
                detailed_msg += "Could not calculate reflectivity for %s\n  %s\n\n" % (xs, sys.exc_value)
                logging.error("Could not calculate reflectivity for %s\n  %s", xs, sys.exc_value)
        return has_errors, detailed_msg

    def calculate_gisans(self, direct_beam):
        has_errors = False
        detailed_msg = ""
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].gisans(direct_beam=direct_beam)
            except:
                raise
                has_errors = True
                detailed_msg += "Could not calculate off-specular reflectivity for %s\n  %s\n\n" % (xs, sys.exc_value)
                logging.error("Could not calculate off-specular reflectivity for %s\n  %s", xs, sys.exc_value)
        if has_errors:
            raise RuntimeError(detailed_msg)

    def calculate_offspec(self, direct_beam=None):
        """
            Loop through the cross-section data sets and update
            the reflectivity.
        """
        has_errors = False
        detailed_msg = ""
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].offspec(direct_beam=direct_beam)
            except:
                has_errors = True
                detailed_msg += "Could not calculate off-specular reflectivity for %s\n  %s\n\n" % (xs, sys.exc_value)
                logging.error("Could not calculate off-specular reflectivity for %s\n  %s", xs, sys.exc_value)
        return has_errors, detailed_msg
            
    def update_configuration(self, configuration):
        """
            Loop through the cross-section data sets and update
            the reflectivity.
        """
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].update_configuration(configuration)
            except:
                logging.error("Could not update configuration for %s\n  %s", xs, sys.exc_value)

    def update_calculated_values(self):
        """
            Loop through the cross-section data sets and update.
        """
        for xs in self.cross_sections:
                self.cross_sections[xs].update_calculated_values()

    def filter_events(self, progress=None):
        """
            Returns a workspace with the selected events
            TODO: make this flexible so we can filter anything we want.

            A cross-section editor dialog should be written to
            add (title, definition) pair list of channels.

            :param function progress: call-back function to track progress
        """
        progress(5, "Preparing...", out_of=100.0)
        nxs = h5py.File(self.file_path, mode='r')
        keys = nxs.keys()
        keys.sort()
        nxs.close()
        # If we have one entry, filter the events to extract channels
        #TODO: for now, just assume it's unpolarized
        self.cross_sections = OrderedDict()
        if len(keys) == 1:
            pass
        else:
            progress_value = 0
            for channel in keys:
                if progress is not None:
                    progress_value += int(100.0/len(keys))
                    progress(progress_value, "Loading %s..." % str(channel), out_of=100.0)
                logging.warning("Loading file %s [%s]", str(self.file_path), str(channel))
                try:
                    nxs_data = LoadEventNexus(str(self.file_path),
                                              #OutputWorkspace="nxs_data",
                                              NXentryName=str(channel))
                except:
                    logging.error("Could not load file %s [%s]", str(self.file_path), str(channel))
                    continue

                name = XSECT_MAPPING.get(channel, 'x')
                cross_section = CrossSectionData(name, self.configuration, entry_name=channel)
                cross_section.collect_info(nxs_data)
                cross_section.process_data(nxs_data)
                self.cross_sections[name] = cross_section
                self.number = cross_section.number
                # Delete workspace
                DeleteWorkspace(nxs_data)

    def load(self, progress=None):
        """
            Load cross-sections from a nexus file.
            :param function progress: call-back function to track progress
        """
        self.filter_events(progress=progress)
        return self.cross_sections

class CrossSectionData(object):
    """
        Data object to hold loaded reflectivity data
    """
    from_event_mode=True
    def __init__(self, name, configuration, entry_name='entry'):
        self.name = name
        self.entry_name = entry_name
        self.measurement_type = 'polarized'
        self.configuration = copy.deepcopy(configuration)
        self.number = 0
        self.q = None
        self._r = None
        self._dr = None
        # Flag to tell us whether we succeeded in using the meta data ROI
        self.use_roi_actual = True
        # Flag to tell us whether we found this data to be a direct beam data set
        self.is_direct_beam = False
        self.tof_range = [0, 0]

    ################## Properties for easy data access ##########################
    # return the size of the data stored in memory for this dataset
    @property
    def nbytes(self): return self.xydata.nbytes+self.xtofdata.nbytes

    @property
    def xdata(self): return self.xydata.mean(axis=0)

    @property
    def ydata(self): return self.xydata.mean(axis=1)

    @property
    def tofdata(self): return self.xtofdata.mean(axis=0)

    # coordinates corresponding to the data items
    @property
    def x(self): return np.arange(self.xydata.shape[1])

    @property
    def y(self): return np.arange(self.xydata.shape[0])

    @property
    def xy(self): return np.meshgrid(self.x, self.y)

    @property
    def tof(self): return (self.tof_edges[:-1]+self.tof_edges[1:])/2.

    @property
    def xtof(self): return np.meshgrid(self.tof, self.x)

    @property
    def r(self):
        if self._r is None:
            return None
        return self._r * self.configuration.scaling_factor

    @property
    def dr(self):
        if self._dr is None:
            return None
        return self._dr * self.configuration.scaling_factor

    @property
    def wavelength(self):
        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        v_n=self.dist_mod_det/self.tof*1e6 #m/s
        return h/m/v_n*1e10 #A

    @property
    def active_area_x(self):
        if self._active_area_x is None:
            return (0, self.xydata.shape[1])
        else:
            return self._active_area_x
    @active_area_x.setter
    def active_area_x(self, value):
        self._active_area_x=value

    @property
    def active_area_y(self):
        if self._active_area_y is None:
            return (0, self.xydata.shape[1])
        else:
            return self._active_area_y

    @active_area_y.setter
    def active_area_y(self, value):
        self._active_area_y=value

    ################## Properties for easy data access ##########################

    def collect_info(self, workspace):
        """
            Extract meta data from DASLogs.

            TODO: get average of values post filtering so that it truly represents the data
        """
        data = workspace.getRun()
        self.origin=(os.path.abspath(data['filename'].value), 'entry')
        self.logs={}
        self.log_minmax={}
        self.log_units={}

        for motor in data.keys():
            if motor in ['proton_charge', 'frequency', 'Veto_pulse']:
                continue
            item = data[motor]
            try:
                self.log_units[motor]=unicode(item.units, encoding='utf8')
                if item.type == 'string':
                    pass
                    #self.logs[motor] = item.value
                    #self.log_minmax[motor] = (item.value, item.value)
                elif item.type == 'number':
                    self.logs[motor] = np.float64(item.value)
                    self.log_minmax[motor] = (np.float64(item.value), np.float64(item.value))
                else:
                    stats = item.getStatistics()
                    self.logs[motor] = np.float64(stats.mean)
                    self.log_minmax[motor] = (np.float64(stats.minimum), np.float64(stats.maximum))
            except:
                logging.error("Error reading DASLogs %s: %s", motor, sys.exc_value)

        self.proton_charge=data['gd_prtn_chrg'].value
        self.total_counts=workspace.getNumberEvents()
        self.total_time=data['duration'].value

        self.experiment=str(data['experiment_identifier'].value)
        self.number=int(workspace.getRunNumber())
        self.merge_warnings=''

        # Retrieve instrument-specific information
        self.configuration.instrument.get_info(workspace, self)

    def process_data(self, workspace):
        run_object = workspace.getRun()

        # Bin events
        if self.configuration.tof_overwrite is not None:
            tof_edges = self.configuration.tof_overwrite
        else:
            tmin, tmax = self.configuration.instrument.get_tof_range(run_object)

            if self.configuration.tof_bin_type == 1: # constant Q
                tof_edges = 1./np.linspace(1./tmin, 1./tmax, self.configuration.tof_bins+1)
            elif self.configuration.tof_bin_type == 2: # constant 1/wavelength
                tof_edges = tmin*(((tmax/tmin)**(1./self.configuration.tof_bins))**np.arange(self.configuration.tof_bins+1))
            else:
                tof_edges = np.linspace(tmin, tmax, self.configuration.tof_bins+1)

        binning_ws = CreateWorkspace(DataX=tof_edges, DataY=np.zeros(len(tof_edges)-1))
        data_rebinned = RebinToWorkspace(WorkspaceToRebin=workspace, WorkspaceToMatch=binning_ws)
        Ixyt = getIxyt(data_rebinned)

        # Create projections for the 2D datasets
        Ixy=Ixyt.sum(axis=2)
        Ixt=Ixyt.sum(axis=1)
        # Store the data
        self.tof_edges=tof_edges
        self.data=Ixyt.astype(float) # 3D dataset
        self.xydata=Ixy.transpose().astype(float) # 2D dataset
        self.xtofdata=Ixt.astype(float) # 2D dataset

        # Determine reduction parameter
        data_info = DataInfo(workspace, self.name, self.configuration)
        self.use_roi_actual = data_info.use_roi_actual
        self.is_direct_beam = data_info.is_direct_beam
        self.tof_range = data_info.tof_range

        self.meta_data_roi_peak = data_info.roi_peak
        self.meta_data_roi_bck = data_info.roi_background

        if not self.configuration.force_peak_roi:
            self.configuration.peak_roi = data_info.peak_range

        if not self.configuration.force_low_res_roi:
            self.configuration.low_res_roi = data_info.low_res_range

        if not self.configuration.force_bck_roi:
            self.configuration.bck_roi = data_info.background

        if self.configuration.set_direct_pixel:
            self.direct_pixel = self.configuration.direct_pixel_overwrite

        if self.configuration.set_direct_angle_offset:
            self.angle_offset = self.configuration.direct_angle_offset_overwrite

        self.scattering_angle = self.configuration.instrument.scattering_angle_from_data(self)

    def update_calculated_values(self):
        """
            Update parameters that are calculated from the configuration
        """
        self.scattering_angle = self.configuration.instrument.scattering_angle_from_data(self)

    def update_configuration(self, configuration):
        """
            Update configuration
        """
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)
            self.update_calculated_values()

    def reflectivity(self, direct_beam=None, configuration=None):
        """
            Compute reflectivity
        """
        self.q = None
        self._r = None
        self._dr = None
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)

        if self.configuration is None:
            return

        # If a direct beam object was passed, use it.
        apply_norm = direct_beam is not None and not self.is_direct_beam
        if not apply_norm:
            direct_beam = CrossSectionData('none', self.configuration, 'none')

        logging.error("%s:%s Reduction with DB: %s [config: %s]",
                      self.number, self.entry_name, direct_beam.number,
                      self.configuration.normalization)
        angle_offset = 0 # Offset from dangle0, in radians
        def _as_ints(a): return [int(a[0]), int(a[1])]
        ws = MagnetismReflectometryReduction(RunNumbers=[str(self.number),],
                                    NormalizationRunNumber=str(direct_beam.number),
                                    SignalPeakPixelRange=_as_ints(self.configuration.peak_roi),
                                    SubtractSignalBackground=True,
                                    SignalBackgroundPixelRange=_as_ints(self.configuration.bck_roi),
                                    ApplyNormalization=apply_norm,
                                    NormPeakPixelRange=_as_ints(direct_beam.configuration.peak_roi),
                                    SubtractNormBackground=True,
                                    NormBackgroundPixelRange=_as_ints(direct_beam.configuration.bck_roi),
                                    CutLowResDataAxis=True,
                                    LowResDataAxisPixelRange=_as_ints(self.configuration.low_res_roi),
                                    CutLowResNormAxis=True,
                                    LowResNormAxisPixelRange=_as_ints(direct_beam.configuration.low_res_roi),
                                    CutTimeAxis=True,
                                    QMin=0.001,
                                    QStep=-0.01,
                                    AngleOffset = angle_offset,
                                    UseWLTimeAxis=False,
                                    TimeAxisStep=self.configuration.tof_bins,
                                    UseSANGLE=not self.configuration.use_dangle,
                                    TimeAxisRange=self.tof_range,
                                    SpecularPixel=self.configuration.peak_position,
                                    ConstantQBinning=self.configuration.use_constant_q,
                                    EntryName=str(self.entry_name))

        self.q = ws.readX(0)[:].copy()
        self._r = ws.readY(0)[:].copy() #* self.configuration.scaling_factor
        self._dr = ws.readE(0)[:].copy() #* self.configuration.scaling_factor
        DeleteWorkspace(ws)

    def offspec(self, direct_beam=None):
        """
            Extract off-specular scattering from 4D dataset (x,y,ToF,I).
            Uses a window in y to filter the 4D data
            and than sums all I values for each ToF and x channel.
            Qz,Qx,kiz,kfz is calculated using the x and ToF positions
            together with the tth-bank and direct pixel values.

            :param CrossSectionData direct_beam: if given, this data will be used to normalize the output
        """
        #TODO: correct for detector sensitivity
        #TODO: Background calculation

        x_pos = self.configuration.peak_position
        x_width = self.configuration.peak_width
        y_pos = self.configuration.low_res_position
        y_width = self.configuration.low_res_width
        scale = 1./self.proton_charge * self.configuration.scaling_factor

        # Get regions in pixels as integers
        reg=map(lambda item: int(round(item)),
                [x_pos-x_width/2., x_pos+x_width/2.+1,
                 y_pos-y_width/2., y_pos+y_width/2.+1])

        rad_per_pixel = self.det_size_x / self.dist_sam_det / self.xydata.shape[1]

        xtth = self.direct_pixel - np.arange(self.data.shape[0])[self.active_area_x[0]:
                                                                 self.active_area_x[1]]
        pix_offset_spec = self.direct_pixel - x_pos
        delta_dangle = self.dangle - self.dangle0
        tth_spec = delta_dangle * np.pi/180. + pix_offset_spec * rad_per_pixel
        af = delta_dangle * np.pi/180. + xtth * rad_per_pixel - tth_spec/2.
        ai = np.ones_like(af) * tth_spec / 2.

        #self._calc_bg()

        v_edges = self.dist_mod_det/self.tof_edges * 1e6 #m/s
        lambda_edges = H_OVER_M_NEUTRON / v_edges * 1e10 #A

        wl = (lambda_edges[:-1] + lambda_edges[1:]) / 2.
        # The resolution for lambda is digital range with equal probability
        # therefore it is the bin size divided by sqrt(12)
        self.d_wavelength = np.abs(lambda_edges[:-1] - lambda_edges[1:]) / np.sqrt(12)
        k = 2. * np.pi / wl

        # calculate reciprocal space, incident and outgoing perpendicular wave vectors
        self.Qz=k[np.newaxis, :]*(np.sin(af)+np.sin(ai))[:, np.newaxis]
        self.Qx=k[np.newaxis, :]*(np.cos(af)-np.cos(ai))[:, np.newaxis]
        self.ki_z=k[np.newaxis, :]*np.sin(ai)[:, np.newaxis]
        self.kf_z=k[np.newaxis, :]*np.sin(af)[:, np.newaxis]

        # calculate ROI intensities and normalize by number of points
        raw_multi_dim=self.data[self.active_area_x[0]:self.active_area_x[1], reg[2]:reg[3], :]
        raw = raw_multi_dim.sum(axis=1)
        d_raw = np.sqrt(raw)

        # normalize data by width in y and multiply scaling factor
        self.intensity = raw/(reg[3]-reg[2]) * scale
        self.d_intensity = d_raw/(reg[3]-reg[2]) * scale
        self.S = self.intensity #self.intensity-self.BG[np.newaxis, :]
        self.dS = self.d_intensity #np.sqrt(self.d_intensity**2+(self.dBG**2)[np.newaxis, :])

        if direct_beam is not None:
            if not direct_beam.configuration.tof_bins == self.configuration.tof_bins:
                logging.error("Trying to normalize with a direct beam data set with different binning")

            norm_raw_multi_dim=direct_beam.data[self.active_area_x[0]:self.active_area_x[1], reg[2]:reg[3], :]
            norm_raw = norm_raw_multi_dim.sum(axis=0).sum(axis=0)
            norm_d_raw = np.sqrt(norm_raw)
            norm_raw /= (reg[3]-reg[2]) * direct_beam.proton_charge * direct_beam.configuration.scaling_factor
            norm_d_raw /= (reg[3]-reg[2]) * direct_beam.proton_charge * direct_beam.configuration.scaling_factor

            idxs=norm_raw>0.
            self.dS[:, idxs]=np.sqrt(
                         (self.dS[:, idxs]/norm_raw[idxs][np.newaxis, :])**2+
                         (self.S[:, idxs]/norm_raw[idxs][np.newaxis, :]**2*norm_d_raw[idxs][np.newaxis, :])**2
                         )
            self.S[:, idxs]/=norm_raw[idxs][np.newaxis, :]
            self.S[:, np.logical_not(idxs)]=0.
            self.dS[:, np.logical_not(idxs)]=0.

    def gisans(self, direct_beam=None):
        """
            Compute GISANS

            :param CrossSectionData direct_beam: if given, this data will be used to normalize the output
        """
        #TODO: Perform sensitivity correction
        x_pos = self.configuration.peak_position
        y_pos = self.configuration.low_res_position
        scale = 1./self.proton_charge * self.configuration.scaling_factor

        rad_per_pixel = self.det_size_x / self.dist_sam_det / self.xydata.shape[1]
        xtth = self.direct_pixel - np.arange(self.data.shape[0])[self.active_area_x[0]:
                                                                 self.active_area_x[1]]
        pix_offset_spec = self.direct_pixel - x_pos
        delta_dangle = self.dangle - self.dangle0
        tth_spec = delta_dangle * np.pi/180. + pix_offset_spec * rad_per_pixel
        af = delta_dangle * np.pi/180. + xtth * rad_per_pixel - tth_spec/2.
        ai = np.ones_like(af) * tth_spec / 2.

        phi=(np.arange(self.data.shape[1])[self.active_area_y[0]:
                                           self.active_area_y[1]]-y_pos)*rad_per_pixel

        v_edges=self.dist_mod_det/self.tof_edges*1e6 #m/s
        lambda_edges = H_OVER_M_NEUTRON/v_edges*1e10 #A
        wl = (lambda_edges[:-1] + lambda_edges[1:]) / 2.
        k = 2.*np.pi / wl

        # calculate ROI intensities and normalize by number of points
        P0 = self.configuration.cut_first_n_points
        PN = len(self.tof) - self.configuration.cut_last_n_points

        # calculate reciprocal space, incident and outgoing perpendicular wave vectors
        Qy=k[np.newaxis, np.newaxis, P0:PN]*(np.sin(phi)*np.cos(af)[:, np.newaxis])[:, :, np.newaxis]
        p_i=k[np.newaxis, np.newaxis, P0:PN]*((0*phi)+np.sin(ai)[:, np.newaxis])[:, :, np.newaxis]
        p_f=k[np.newaxis, np.newaxis, P0:PN]*((0*phi)+np.sin(af)[:, np.newaxis])[:, :, np.newaxis]
        Qz=p_i+p_f

        raw=self.data[self.active_area_x[0]:self.active_area_x[1],
                        self.active_area_y[0]:self.active_area_y[1],
                        P0:PN]

        intensity = scale * np.array(raw)
        d_intensity = scale * np.sqrt(raw)

        if direct_beam is not None:
            if not direct_beam.configuration.tof_bins == self.configuration.tof_bins:
                logging.error("Trying to normalize with a direct beam data set with different binning")

            norm_raw_multi_dim=direct_beam.data[self.active_area_x[0]:self.active_area_x[1],
                                                self.active_area_y[0]:self.active_area_y[1], P0:PN]
            norm_raw = norm_raw_multi_dim.sum(axis=0).sum(axis=0)
            norm_d_raw = np.sqrt(norm_raw)
            
            surface = (self.active_area_x[1]-self.active_area_x[0]) * (self.active_area_y[1]-self.active_area_y[0])
            
            norm_raw /= surface * direct_beam.proton_charge * direct_beam.configuration.scaling_factor
            norm_d_raw /= surface * direct_beam.proton_charge * direct_beam.configuration.scaling_factor

            idxs=norm_raw>0.
            d_intensity[:, :, idxs]=np.sqrt(
                         (d_intensity[:, :, idxs]/norm_raw[idxs][np.newaxis, np.newaxis, :])**2+
                         (intensity[:, :, idxs]/norm_raw[idxs][np.newaxis, np.newaxis, :]**2*norm_d_raw[idxs][np.newaxis, np.newaxis, :])**2
                         )
            intensity[:, :, idxs]/=norm_raw[idxs][np.newaxis, np.newaxis, :]
            intensity[:, :, np.logical_not(idxs)]=0.
            d_intensity[:, :, np.logical_not(idxs)]=0.

        # Create grid
        # bins=(self.options['gisans_gridy'], self.options['gisans_gridz']),
        #TODO: allow binning as application parameter
        self.SGrid, qy, qz = np.histogram2d(Qy.flatten(), Qz.flatten(),
                                            bins=(50, 50),
                                            weights=intensity.flatten())
        npoints, _, _ = np.histogram2d(Qy.flatten(), Qz.flatten(),
                                       bins=(50, 50))
        self.SGrid[npoints>0]/=npoints[npoints>0]
        self.SGrid=self.SGrid.transpose()
        qy=(qy[:-1]+qy[1:])/2.
        qz=(qz[:-1]+qz[1:])/2.
        self.QyGrid, self.QzGrid = np.meshgrid(qy, qz)

class NexusMetaData(object):
    """
        Class used to hold meta-data read before loading the neutron events
    """
    mid_q = 0
    is_direct_beam = False

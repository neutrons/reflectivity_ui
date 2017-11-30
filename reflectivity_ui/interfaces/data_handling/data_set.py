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

    def set_parameter(self, param, value):
        """
            Loop through the cross-section data sets and update
            a parameter.
        """
        try:
            for xs in self.cross_sections:
                if hasattr(self.cross_sections[xs].configuration, param):
                    setattr(self.cross_sections[xs].configuration, param, value)
        except:
            logging.error("Could not set parameter %s %s\n  %s", param, value, sys.exc_value)

    def calculate_reflectivity(self, direct_beam=None, configuration=None):
        """
            Loop through the cross-section data sets and update
            the reflectivity.
        """
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].reflectivity(direct_beam=direct_beam, configuration=configuration)
            except:
                logging.error("Could not calculate reflectivity for %s\n  %s", xs, sys.exc_value)

    def filter_events(self):
        """
            Returns a workspace with the selected events
            TODO: make this flexible so we can filter anything we want.

            A cross-section editor dialog should be written to
            add (title, definition) pair list of channels.
        """
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
            for channel in keys:
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

    def load(self):
        """
            Load cross-sections from a nexus file.
        """
        self.filter_events()
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
        self.r = None
        self.dr = None
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
        logging.error(str(data_info.roi_peak))

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

        self.scattering_angle = self.configuration.instrument.scattering_angle(workspace,
                                                                               self.configuration.peak_position,
                                                                               self.direct_pixel,
                                                                               self.angle_offset)

        #self.reflectivity()

    def reflectivity(self, direct_beam=None, configuration=None):
        """
            Compute reflectivity
        """
        if configuration is not None:
            self.configuration = copy.deepcopy(configuration)

        if self.configuration is None:
            return

        # If a direct beam object was passed, use it.
        logging.error("Config DB: %s", self.configuration.normalization)
        apply_norm = direct_beam is not None
        if not apply_norm:
            direct_beam = CrossSectionData('none', self.configuration, 'none')

        logging.error("Reduction with DB: %s", direct_beam.number)
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
        self.r = ws.readY(0)[:].copy()
        self.dr = ws.readE(0)[:].copy()
        DeleteWorkspace(ws)

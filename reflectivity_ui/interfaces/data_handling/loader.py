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

# Import mantid according to the application configuration
from ..configuration import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *

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
        self.configuration = configuration
        self.cross_sections = {}

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
                cross_section = CrossSectionData(name, self.configuration)
                cross_section.collect_info(nxs_data)
                cross_section.process_data(nxs_data)
                self.cross_sections[name] = cross_section

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
    def __init__(self, name, configuration):
        self.name = name
        self.measurement_type = 'polarized'
        self.configuration = configuration

    ################## Properties for easy data access ##########################
    # return the size of the data stored in memory for this dataset
    @property
    def nbytes(self): return len(self._data_zipped)+
                              self.xydata.nbytes+self.xtofdata.nbytes
    @property
    def rawbytes(self): return self.data.nbytes+self.xydata.nbytes+self.xtofdata.nbytes

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
    def lamda(self):
        h = 6.626e-34  # m^2 kg s^-1
        m = 1.675e-27  # kg
        v_n=self.dist_mod_det/self.tof*1e6 #m/s
        lamda_n=h/m/v_n*1e10 #A
        return lamda_n

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
        self.number=int(data['run_number'].value)
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

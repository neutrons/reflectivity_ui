"""
    Loader for event nexus files.
    Uses Mantid Framework
"""
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

class NexusData(object):
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

                run_object = nxs_data.getRun()

                # Bin events
                tmin, tmax = self.configuration.instrument.get_tof_range(run_object)

                if self.configuration.tof_bin_type == 1: # constant Q
                    tof_edges = 1./np.linspace(1./tmin, 1./tmax, self.configuration.tof_bins+1)
                elif self.configuration.tof_bin_type == 2: # constant 1/wavelength
                    tof_edges = tmin*(((tmax/tmin)**(1./self.configuration.tof_bins))**np.arange(self.configuration.tof_bins+1))
                else:
                    tof_edges = np.linspace(tmin, tmax, self.configuration.tof_bins+1)
        
                name = XSECT_MAPPING.get(channel, 'x')
                self.cross_sections[name] = CrossSectionData(name, nxs_data)
        
    def load(self):
        """
            Load cross-sections from a nexus file.
        """
        self.filter_events()
        return self.cross_sections

class CrossSectionData(object):
    """
    """
    def __init__(self, name='x', workspace=None):
        self.name = name
        self.measurement_type = 'polarized'
        self.workspace = workspace
        if workspace is not None:
            self.collect_info()
    
    def collect_info(self):
        """
            Extract meta data from DASLogs.

            TODO: get average of values post filtering so that it truly represents the data
        """
        data = self.workspace.getRun()
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
    
        self.lambda_center=data['LambdaRequest'].value[0]
        self.dangle=data['DANGLE'].value[0]
        if 'DANGLE0' in data: # compatibility for ancient file format
            self.dangle0=data['DANGLE0'].value[0]
            self.dpix=data['DIRPIX'].value[0]
            self.slit1_width=data['S1HWidth'].value[0]
            self.slit2_width=data['S2HWidth'].value[0]
            self.slit3_width=data['S3HWidth'].value[0]
        else:
            self.slit1_width=data['RSlit1'].value[0]-data['LSlit1'].value[0]
            self.slit2_width=data['RSlit2'].value[0]-data['LSlit2'].value[0]
            self.slit3_width=data['RSlit3'].value[0]-data['LSlit3'].value[0]
    
        #TODO: these don't exist in the DASLogs
        #self.slit1_dist=-data['instrument/aperture1/distance'].value[0]*1000.
        #self.slit2_dist=-data['instrument/aperture2/distance'].value[0]*1000.
        #self.slit3_dist=-data['instrument/aperture3/distance'].value[0]*1000.
    
        self.sangle=data['SANGLE'].value[0]
        self.proton_charge=data['gd_prtn_chrg'].value
        self.total_counts=self.workspace.getNumberEvents()
        self.total_time=data['duration'].value
    
        self.dist_sam_det=data['SampleDetDis'].value[0]*1e-3
        self.dist_mod_det=data['ModeratorSamDis'].value[0]*1e-3+self.dist_sam_det
        self.dist_mod_mon=data['ModeratorSamDis'].value[0]*1e-3-2.75
        
        # Get these from instrument
        self.det_size_x = int(self.workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0]) #304
        self.det_size_y = int(self.workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0]) #256
    
        self.experiment=str(data['experiment_identifier'].value)
        self.number=int(data['run_number'].value)
        self.merge_warnings=''
    
        # The following active area used to be taken from instrument.DETECTOR_REGION
        self.active_area_x = (8, 295)
        self.active_area_y = (8, 246)
    

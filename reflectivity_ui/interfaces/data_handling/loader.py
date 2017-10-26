"""
    Loader for event nexus files.
    Uses Mantid Framework
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
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
                self.cross_sections[name] = CrossSectionData(name)
        
    def load(self):
        """
            Load cross-sections from a nexus file.
        """
        self.filter_events()
        return self.cross_sections

class CrossSectionData(object):
    """
    """
    def __init__(self, name='x'):
        self.name = name
        
        
        
        
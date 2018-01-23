"""
    Methods used to process data, usually calling Mantid
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long, multiple-statements
from __future__ import absolute_import, division, print_function
import sys
import logging
import h5py

# Import mantid according to the application configuration
from . import ApplicationConfiguration
from ase.cluster.factory import cross
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *

from .instrument import Instrument
from .data_set import NexusMetaData


def stitch_reflectivity(reduction_list, xs=None, normalize_to_unity=True):
    """
        Stitch and normalize data sets
        
        :param string xs: name of the cross-section to use
        :param bool normalize_to_unity: if True, the specular ridge will be normalized to 1
    """
    if len(reduction_list) == 0:
        return []

    # Select the cross-section we will use to determine the scaling factors
    if xs is None:
        xs = reduction_list[0].cross_sections.keys()[0]

    # First, determine the overall scaling factor as needed
    scaling_factor = 1.0
    if normalize_to_unity:
        idx_list = reduction_list[0].cross_sections[xs].q < reduction_list[0].cross_sections[xs].configuration.total_reflectivity_q_cutoff
        total = 0
        weights = 0
        for i in range(len(reduction_list[0].cross_sections[xs]._r)):
            if idx_list[i]:
                w = 1.0 / float(reduction_list[0].cross_sections[xs]._dr[i])**2
                total += w * float(reduction_list[0].cross_sections[xs]._r[i])
                weights += w
        if weights > 0:
            scaling_factor = weights / total
        reduction_list[0].set_parameter("scaling_factor", scaling_factor)

    # Stitch the data sets together
    _previous_ws = None
    running_scale = scaling_factor
    scaling_factors = [running_scale]

    for i in range(len(reduction_list)):
        ws = CreateWorkspace(DataX=reduction_list[i].cross_sections[xs].q,
                             DataY=reduction_list[i].cross_sections[xs]._r,
                             DataE=reduction_list[i].cross_sections[xs]._dr)
        ws = ConvertToHistogram(ws)
        if _previous_ws is not None:
            _, scale = Stitch1D(_previous_ws, ws)
            running_scale *= scale
            scaling_factors.append(running_scale)
            reduction_list[i].set_parameter("scaling_factor", running_scale)
        _previous_ws = CloneWorkspace(ws)

    return scaling_factors

def extract_meta_data(file_path=None, cross_section_data=None, configuration=None):
    """
        Get mid Q-value from meta data
        :param str file_path: name of the file to read
    """
    meta_data = NexusMetaData()

    if cross_section_data is not None:
        meta_data.mid_q = cross_section_data.configuration.instrument.mid_q_value_from_data(cross_section_data)
        meta_data.is_direct_beam = cross_section_data.is_direct_beam
        return meta_data
    elif file_path is None:
        raise RuntimeError("Either a file path or a data object must be supplied")

    nxs = h5py.File(file_path, mode='r')
    keys = nxs.keys()
    keys.sort()
    nxs.close()

    if len(keys) == 0:
        logging.error("No entry in data file %s", file_path)
        return meta_data

    try:
        ws = LoadEventNexus(str(file_path),
                            MetaDataOnly=True,
                            NXentryName=str(keys[0]))
        meta_data.mid_q = Instrument.mid_q_value(ws)
        if configuration is not None:
            meta_data.is_direct_beam = configuration.instrument.check_direct_beam(ws)
    except:
        logging.error(sys.exc_value)
        logging.error("Could not load file %s [%s]", file_path, keys[0])

    return meta_data

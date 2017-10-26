#-*- coding: utf-8 -*-
"""
    This code is used to bypass the QuickNXS loader and use Mantid instead.
        
    The code should be hooked up in the NXSData class by adding the following
    at the top of the _read_file() method:
    
    def _read_file(self, filename):
      if filename.endswith('nxs.h5'):
          import mantid_utils
          return mantid_utils.read_file(self, filename, MRDataset)
        
"""
from __future__ import division
import sys
import os
import time
import h5py
import logging
import numpy as np
sys.path.insert(0,'/opt/Mantid/bin')
from mantid.simpleapi import *

### Parameters needed for some calculations.
H_OVER_M_NEUTRON=3.956034e-7 # h/m_n [m²/s]
ANALYZER_IN=(0., 100.) # position and maximum deviation of analyzer in it's working position
POLARIZER_IN=(-348., 50.) # position and maximum deviation of polarizer in it's working position
SUPERMIRROR_IN=(19.125, 10.) # position and maximum deviation of the supermirror translation
POLY_CORR_PARAMS=[-4.74152261e-05,-4.62469580e-05, 1.25995446e-02, 2.13654008e-02,
                  1.02334517e+01] # parameters used in polynomial detector sensitivity correction
DETECTOR_SENSITIVITY={}
# measurement type mapping of states
MAPPING_12FULL=(
                 (u'++ (0V)', u'entry-off_off_Ezero'),
                 (u'-- (0V)', u'entry-on_on_Ezero'),
                 (u'+- (0V)', u'entry-off_on_Ezero'),
                 (u'-+ (0V)', u'entry-on_off_Ezero'),
                 (u'++ (+V)', u'entry-off_off_Eplus'),
                 (u'-- (+V)', u'entry-on_on_Eplus'),
                 (u'+- (+V)', u'entry-off_on_Eplus'),
                 (u'-+ (+V)', u'entry-on_off_Eplus'),
                 (u'++ (-V)', u'entry-off_off_Eminus'),
                 (u'-- (-V)', u'entry-on_on_Eminus'),
                 (u'+- (-V)', u'entry-off_on_Eminus'),
                 (u'-+ (-V)', u'entry-on_off_Eminus'),
                 )
MAPPING_12HALF=(
                 (u'+ (0V)', u'entry-off_off_Ezero'),
                 (u'- (0V)', u'entry-on_off_Ezero'),
                 (u'+ (+V)', u'entry-off_off_Eplus'),
                 (u'- (+V)', u'entry-on_off_Eplus'),
                 (u'+ (-V)', u'entry-off_off_Eminus'),
                 (u'- (-V)', u'entry-on_off_Eminus'),
                 )
MAPPING_FULLPOL=(
                 (u'++', u'entry-Off_Off'),
                 (u'--', u'entry-On_On'),
                 (u'+-', u'entry-Off_On'),
                 (u'-+', u'entry-On_Off'),
                 )
MAPPING_HALFPOL=(
                 (u'+', u'entry-Off_Off'),
                 (u'-', u'entry-On_Off'),
                 )
MAPPING_UNPOL=(
               (u'x', u'entry-Off_Off'),
               )
MAPPING_EFIELD=(
                (u'0V', u'entry-Off_Off'),
                (u'+V', u'entry-On_Off'),
                (u'-V', u'entry-Off_On'),
                )

try:
    from .ipython_tools import NiceDict
except:
    def NiceDict():
        return {}


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
    _y_error_axis = np.zeros((sz_x_axis, sz_y_axis, nbr_tof-1))

    for x in range(sz_x_axis):
        for y in range(sz_y_axis):
            _index = int(sz_y_axis*x+y)
            _tmp_data = nxs_data.readY(_index)[:]
            _y_axis[x,y,:] = _tmp_data

    return _y_axis

# Method of class MRDataset
def collect_info(self, nxs_data):
    """
        Extract meta data from DASLogs.
        @param self: MRDataset object
        @param nxs_data: Mantid workspace
    """
    data = nxs_data.getRun()
    self.origin=(os.path.abspath(data['filename'].value), 'entry')
    self.logs=NiceDict()
    self.log_minmax=NiceDict()
    self.log_units=NiceDict()

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
    self.total_counts=nxs_data.getNumberEvents()
    self.total_time=data['duration'].value

    self.dist_sam_det=data['SampleDetDis'].value[0]*1e-3
    self.dist_mod_det=data['ModeratorSamDis'].value[0]*1e-3+self.dist_sam_det
    self.dist_mod_mon=data['ModeratorSamDis'].value[0]*1e-3-2.75
    
    # Get these from instrument
    self.det_size_x = int(nxs_data.getInstrument().getNumberParameter("number-of-x-pixels")[0]) #304
    self.det_size_y = int(nxs_data.getInstrument().getNumberParameter("number-of-y-pixels")[0]) #256

    self.experiment=str(data['experiment_identifier'].value[0])
    self.number=int(data['run_number'].value[0])
    self.merge_warnings=''

    # The following active area used to be taken from instrument.DETECTOR_REGION
    self.active_area_x = (8, 295)
    self.active_area_y = (8, 246)


# Method of class MRDataset
def from_nxs_h5(cls, data, read_options,
                callback=None, callback_offset=0., callback_scaling=1.,
                tof_overwrite=None):
    """
        Load data from a Nexus file containing event information.
        Creates 3D histogram with ither linear or 1/t spaced 
        time of flight channels. The result has the same format as
        from the read_file function.
    """
    output=cls()
    output.read_options=read_options
    output.from_event_mode=True
    bin_type=read_options['bin_type']
    bins=read_options['bins']
    try:
        collect_info(output, data)
    except KeyError:
        logging.warn('Error while collecting metadata:\n\n'+traceback.format_exc())

    run_object = data.getRun()
    if tof_overwrite is None:
        lcenter=run_object['LambdaRequest'].value[0]
        # ToF region for this specific central wavelength
        tmin=output.dist_mod_det/H_OVER_M_NEUTRON*(lcenter-1.6)*1e-4
        tmax=output.dist_mod_det/H_OVER_M_NEUTRON*(lcenter+1.6)*1e-4

        if bin_type==0: # constant ??
            tof_edges=np.linspace(tmin, tmax, bins+1)
        elif bin_type==1: # constant ?Q
            tof_edges=1./np.linspace(1./tmin, 1./tmax, bins+1)
        elif bin_type==2: # constant ??/?
            tof_edges=tmin*(((tmax/tmin)**(1./bins))**np.arange(bins+1))
        else:
            raise ValueError('Unknown bin type %i' % bin_type)
    else:
        tof_edges=tof_overwrite

    #TODO: Not sure what the split bins is supposed to mean
    split_bins=read_options['event_split_bins']
    split_index=read_options['event_split_index']
    
    binning_ws = CreateWorkspace(DataX=tof_edges, DataY=np.zeros(len(tof_edges)-1))
    data_rebinned = RebinToWorkspace(WorkspaceToRebin=data, WorkspaceToMatch=binning_ws)
    Ixyt = getIxyt(data_rebinned)

    # Create projections for the 2D datasets
    Ixy=Ixyt.sum(axis=2)
    Ixt=Ixyt.sum(axis=1)
    # Store the data
    output.tof_edges=tof_edges
    output.data=Ixyt.astype(float) # 3D dataset
    output.xydata=Ixy.transpose().astype(float) # 2D dataset
    output.xtofdata=Ixt.astype(float) # 2D dataset
    return output


# Method of class NXSData
def read_file(self, filename, data_cls):
    """
        Load data from a Nexus file.
    """
    start=time.time()
    if self._options['callback']:
        self._options['callback'](0.)
    try:
        nxs_data = LoadEventNexus(Filename=str(filename), OutputWorkspace="nxs_data", MetaDataOnly=True)
        nxs=h5py.File(filename, mode='r')
    except IOError:
        logging.error('Could not read nxs file %s'%filename, exc_info=True)
        return False
        
    # Analyze channels
    channels=nxs.keys()
    try:
        # Here we should filter by polarization and electric field
        max_counts = max([nxs[channel][u'total_counts'].value[0] for channel in channels])
    except KeyError:
        logging.error('total_counts not defined in channels')
        return False
    for channel in list(channels):
        if nxs[channel][u'total_counts'].value[0]<(self.COUNT_THREASHOLD*max_counts):
            channels.remove(channel)
    if len(channels)==0:
        logging.error('No valid channels in file')
        return False

    run_object = nxs_data.getRun()
    ana = run_object['AnalyzerLift'].value[0]
    pol = run_object['PolLift'].value[0]
    smpt = run_object['SMPolTrans'].value[0]

    # Select the type of measurement that has been used
    #TODO: This should be replaced by analyzing the logs
    if abs(ana-ANALYZER_IN[0])<ANALYZER_IN[1]: # is analyzer is in position
        if channels[0] in [m[1] for m in MAPPING_12FULL]:
            self.measurement_type='Polarization Analysis w/E-Field'
            mapping=list(MAPPING_12FULL)
        else:
            self.measurement_type='Polarization Analysis'
            mapping=list(MAPPING_FULLPOL)
    elif abs(pol-POLARIZER_IN[0])<POLARIZER_IN[1] or \
        abs(smpt-SUPERMIRROR_IN[0])<SUPERMIRROR_IN[1]: # is bender or supermirror polarizer is in position
        if channels[0] in [m[1] for m in MAPPING_12HALF]:
            self.measurement_type='Polarized w/E-Field'
            mapping=list(MAPPING_12HALF)
        else:
            self.measurement_type='Polarized'
            mapping=list(MAPPING_HALFPOL)
    elif 'DASlogs' in nxs[channels[0]] and \
        run_object['SP_HV_Minus'].getStatistics().standard_deviation > 0:
        #nxs[channels[0]]['DASlogs'].get('SP_HV_Minus') is not None and \
        #channels!=[u'entry-Off_Off']: # is E-field cart connected and not only 0V measured
        self.measurement_type='Electric Field'
        mapping=list(MAPPING_EFIELD)
    elif len(channels)==1:
        self.measurement_type='Unpolarized'
        mapping=list(MAPPING_UNPOL)
    else:
        self.measurement_type='Unknown'
        mapping=[]
    # check that all channels have a mapping entry
    for channel in channels:
        if not channel in [m[1] for m in mapping]:
            mapping.append((channel.lstrip('entry-'), channel))

    progress=0.1
    if self._options['callback']:
        self._options['callback'](progress)
    self._read_times.append(time.time()-start)
    i=1
    empty_channels=[]
    
    for dest, channel in mapping:
        if channel not in channels:
            continue
        raw_data=nxs[channel]
        #TODO: The filtering might be done here
        nxs_data = LoadEventNexus(str(filename), OutputWorkspace="nxs_data", NXentryName=str(channel))
      
        data = from_nxs_h5(data_cls, nxs_data, self._options,
                           callback=self._options['callback'],
                           callback_offset=progress,
                           callback_scaling=0.9/len(channels),
                           tof_overwrite=self._options['event_tof_overwrite'])

        # Update the origin data member with the cross-section channel
        data.origin = (data.origin[0], channel)
        if data is None:
            # no data in channel, don't add it
            empty_channels.append(dest)
            continue

        self._channel_data.append(data)
        self._channel_names.append(dest)
        self._channel_origin.append(channel)
        progress=0.1+0.9*float(i)/len(channels)
        if self._options['callback']:
            self._options['callback'](progress)
        i+=1
        self._read_times.append(time.time()-self._read_times[-1]-start)

    nxs.close()
    if empty_channels:
        logging.warn('No counts for state %s'%(','.join(empty_channels)))
    return True


"""
    Loader for event nexus files.
    Uses Mantid Framework
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long, multiple-statements, bare-except, wrong-import-order, too-many-locals, too-few-public-methods
from __future__ import absolute_import, division, print_function
import sys
import logging
from collections import OrderedDict
import copy
import math
import numpy as np

# Import mantid according to the application configuration
from . import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
import mantid.simpleapi as api
# Set Mantid logging level to warnings
ConfigService.setConsoleLogLevel(4)

from .data_info import DataInfo
from . import off_specular

### Parameters needed for some calculations.

ANALYZER_IN = (0., 100.) # position and maximum deviation of analyzer in it's working position
POLARIZER_IN = (-348., 50.) # position and maximum deviation of polarizer in it's working position
SUPERMIRROR_IN = (19.125, 10.) # position and maximum deviation of the supermirror translation
POLY_CORR_PARAMS = [-4.74152261e-05, -4.62469580e-05, 1.25995446e-02, 2.13654008e-02,
                    1.02334517e+01] # parameters used in polynomial detector sensitivity correction

XSECT_MAPPING = {u'entry-Off_Off': u'++',
                 u'entry-On_On': u'--',
                 u'entry-Off_On': u'+-',
                 u'entry-On_Off': u'-+',
                }

H_OVER_M_NEUTRON = 3.956034e-7 # h/m_n [m^2/s]

# Number of events under which we throw away a workspace
#TODO: This should be a parameter
N_EVENTS_CUTOFF = 100

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
            _y_axis[x, y, :] = _tmp_data

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
        if has_errors:
            raise RuntimeError(detailed_msg)

    def calculate_gisans(self, direct_beam):
        has_errors = False
        detailed_msg = ""
        for xs in self.cross_sections:
            try:
                self.cross_sections[xs].gisans(direct_beam=direct_beam)
            except:
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
        if has_errors:
            raise RuntimeError(detailed_msg)

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

    def map_cross_section(self, ws):
        """
            Return the proper cross-section label for the provided workspace.
            For legacy data, we just map Off to + and On to -.
            #TODO: For new data, use the log entries to determine whether Off is + or -.
            :param workspace ws: workspace to inspect
        """
        entry = ws.getRun().getProperty("cross_section_id").value
        return XSECT_MAPPING.get('entry-%s' % entry, 'x')

    def load(self, progress=None):
        """
            Load cross-sections from a nexus file.
            :param function progress: call-back function to track progress
        """
        if progress is not None:
            progress(5, "Filtering data...", out_of=100.0)

        try:
            xs_list = self.configuration.instrument.load_data(self.file_path)
            logging.info("%s loaded: %s xs", self.file_path, len(xs_list))
        except:
            logging.error("Could not load file %s\n  %s", str(self.file_path), sys.exc_value)
            return

        self.cross_sections = OrderedDict()
        progress_value = 0
        # Keep track of cross-section with max counts so we can use it to
        # select peak regions
        _max_counts = 0
        _max_xs = None
        for ws in xs_list:
            # Get the unique name for the cross-section, determined by the filtering
            channel = ws.getRun().getProperty("cross_section_id").value
            if progress is not None:
                progress_value += int(100.0/len(xs_list))
                progress(progress_value, "Loading %s..." % str(channel), out_of=100.0)

            # Get rid of emty workspaces
            logging.info("Loading %s: %s events", str(channel), ws.getNumberEvents())
            if ws.getNumberEvents() < N_EVENTS_CUTOFF:
                logging.warn("Too few events for %s: %s", channel, ws.getNumberEvents())
                continue

            name = self.map_cross_section(ws)
            cross_section = CrossSectionData(name, self.configuration, entry_name=channel)
            cross_section.collect_info(ws)
            cross_section.process_data(ws)
            self.cross_sections[name] = cross_section
            self.number = cross_section.number
            if cross_section.total_counts > _max_counts:
                _max_counts = cross_section.total_counts
                _max_xs = name
        # Push the configuration (reduction options and peak regions) from the
        # cross-section with the most data to all other cross-sections.
        if _max_xs is not None:
            self.cross_sections[_max_xs].get_reduction_parameters()
            self.update_configuration(self.cross_sections[_max_xs].configuration)

        progress(100, "Complete", out_of=100.0)

        return self.cross_sections


class CrossSectionData(object):
    """
        Data object to hold loaded reflectivity data
    """
    # Instrument specific attributes filled out by the instrument object
    dist_sam_det = 0
    dist_mod_det = 0
    dist_mod_mon = 0
    dist_mod_mon = 0
    dist_mod_mon = 0
    dangle = 0
    dangle0 = 0
    det_size_x = 0.0007
    det_size_y = 0.0007

    def __init__(self, name, configuration, entry_name='entry'):
        self.name = name
        self.entry_name = entry_name
        self.measurement_type = 'polarized'
        self.configuration = copy.deepcopy(configuration)
        self.number = 0
        self.q = None
        self._r = None
        self._dr = None
        self._event_workspace = None
        # Flag to tell us whether we succeeded in using the meta data ROI
        self.use_roi_actual = True
        # Flag to tell us whether we found this data to be a direct beam data set
        self.is_direct_beam = False
        self.tof_range = [0, 0]
        self._active_area_x = None
        self._active_area_y = None
        self.logs = {}
        self.log_minmax = {}
        self.log_units = {}
        self.proton_charge = 0
        self.total_counts = 0
        self.total_time = 0

        self.experiment = ''
        self.number = 0
        self.merge_warnings = ''
        self.tof_edges = None
        self.data = None
        self.xydata = None
        self.xtofdata = None
        self.meta_data_roi_peak = None
        self.meta_data_roi_bck = None
        self.direct_pixel = 0
        self.angle_offset = 0
        self.scattering_angle = 0
        self._reflectivity_workspace = None

        # Offset data
        self.off_spec = None

        # GISANS
        #TODO: refactor this
        self.SGrid = None
        self.QyGrid = None
        self.QzGrid = None

    ################## Properties for easy data access ##########################
    # return the size of the data stored in memory for this dataset
    #pylint: disable=missing-docstring
    @property
    def reflectivity_workspace(self):
        if str(self._reflectivity_workspace) in api.mtd:
            return api.mtd[self._reflectivity_workspace]
        return None

    @property
    def nbytes(self): return self.data.nbytes

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
        v_n = self.dist_mod_det/self.tof*1e6 #m/s
        return h/m/v_n*1e10 #A

    @property
    def active_area_x(self):
        if self._active_area_x is None:
            return (0, self.xydata.shape[1])
        else:
            return self._active_area_x
    @active_area_x.setter
    def active_area_x(self, value):
        self._active_area_x = value

    @property
    def active_area_y(self):
        if self._active_area_y is None:
            return (0, self.xydata.shape[1])
        else:
            return self._active_area_y

    @active_area_y.setter
    def active_area_y(self, value):
        self._active_area_y = value

    #pylint: enable=missing-docstring
    ################## Properties for easy data access ##########################

    def collect_info(self, workspace):
        """
            Extract meta data from DASLogs.

            TODO: get average of values post filtering so that it truly represents the data
        """
        data = workspace.getRun()
        #self.origin=(os.path.abspath(data['filename'].value), 'entry')
        self.logs = {}
        self.log_minmax = {}
        self.log_units = {}

        for motor in data.keys():
            if motor in ['proton_charge', 'frequency', 'Veto_pulse']:
                continue
            item = data[motor]
            try:
                self.log_units[motor] = unicode(item.units, encoding='utf8')
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

        self.proton_charge = data['gd_prtn_chrg'].value
        self.total_counts = workspace.getNumberEvents()
        self.total_time = data['duration'].value

        self.experiment = str(data['experiment_identifier'].value)
        self.number = int(workspace.getRunNumber())
        self.merge_warnings = ''

        # Retrieve instrument-specific information
        self.configuration.instrument.get_info(workspace, self)

    def process_data(self, workspace):
        self._event_workspace = str(workspace)

        # Get initial reduction parameters. They may be overwritten later.
        self.get_reduction_parameters(workspace)

        # Bin events
        if self.configuration.tof_overwrite is not None:
            tof_edges = self.configuration.tof_overwrite
        else:
            if self.configuration.tof_bin_type == 1: # constant Q
                tof_edges = 1./np.linspace(1./self.tof_range[0], 1./self.tof_range[1], self.configuration.tof_bins+1)
            elif self.configuration.tof_bin_type == 2: # constant 1/wavelength
                tof_edges = self.tof_range[0]*(((self.tof_range[1]/self.tof_range[0])**(1./self.configuration.tof_bins))**np.arange(self.configuration.tof_bins+1))
            else:
                tof_edges = np.linspace(self.tof_range[0], self.tof_range[1], self.configuration.tof_bins+1)

        binning_ws = api.CreateWorkspace(DataX=tof_edges, DataY=np.zeros(len(tof_edges)-1))
        data_rebinned = api.RebinToWorkspace(WorkspaceToRebin=workspace, WorkspaceToMatch=binning_ws)
        Ixyt = getIxyt(data_rebinned)

        # Create projections for the 2D datasets
        Ixy = Ixyt.sum(axis=2)
        Ixt = Ixyt.sum(axis=1)
        # Store the data
        self.tof_edges = tof_edges
        self.data = Ixyt.astype(float) # 3D dataset
        self.xydata = Ixy.transpose().astype(float) # 2D dataset
        self.xtofdata = Ixt.astype(float) # 2D dataset

        if self.configuration.set_direct_pixel:
            self.direct_pixel = self.configuration.direct_pixel_overwrite

        if self.configuration.set_direct_angle_offset:
            self.angle_offset = self.configuration.direct_angle_offset_overwrite

        self.scattering_angle = self.configuration.instrument.scattering_angle_from_data(self)

    def get_reduction_parameters(self, workspace=None):
        """
            Determine reduction parameter
            :param workspace: mantid workspace
        """
        if workspace is None:
            if self._event_workspace is None:
                return
            else:
                workspace = api.mtd[self._event_workspace]

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

    def get_counts_vs_TOF(self):
        """
            Used for normalization, returns ROI counts vs TOF.
        """
        # Calculate ROI intensities and normalize by number of points
        raw_data = self.data[self.configuration.peak_roi[0]:self.configuration.peak_roi[1],
                             self.configuration.low_res_roi[0]:self.configuration.low_res_roi[1], :]
        size_roi = float((self.configuration.low_res_roi[1] - self.configuration.low_res_roi[0]) * (self.configuration.peak_roi[1] - self.configuration.peak_roi[0]))
        summed_raw = raw_data.sum(axis=0).sum(axis=0)

        # Remove the background
        bck = self.get_background_vs_TOF()

        return (summed_raw/math.fabs(size_roi) - bck)/self.proton_charge

    def get_background_vs_TOF(self):
        """
            Returns the background counts vs TOF
        """
        # Find the background pixels to use, excluding the peak if there's an overlap.
        dims = self.data.shape
        indices = [(i >= self.configuration.bck_roi[0] and i < self.configuration.bck_roi[1]) \
                   and not (i >= self.configuration.peak_roi[0] and i < self.configuration.peak_roi[1]) for i in range(dims[0])]
        indices = np.asarray(indices)
        n_bins = len(indices[indices == True])
        raw_bck = self.data[indices,
                            self.configuration.low_res_roi[0]:self.configuration.low_res_roi[1], :]
        summed_bck = raw_bck.sum(axis=0).sum(axis=0)
        size_bck = float(n_bins * (self.configuration.low_res_roi[1]-self.configuration.low_res_roi[0]))

        return summed_bck/math.fabs(size_bck)

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
        apply_norm = direct_beam is not None # and not self.is_direct_beam
        if not apply_norm:
            direct_beam = CrossSectionData('none', self.configuration, 'none')

        logging.error("%s:%s Reduction with DB: %s [config: %s]",
                      self.number, self.entry_name, direct_beam.number,
                      self.configuration.normalization)
        angle_offset = 0 # Offset from dangle0, in radians
        def _as_ints(a): return [int(round(a[0])), int(round(a[1]))]
        output_ws = "r%s_%s" % (self.number, str(self.entry_name))

        ws_norm = None
        if apply_norm and direct_beam._event_workspace is not None:
            ws_norm = api.CloneWorkspace(InputWorkspace=direct_beam._event_workspace)

        ws = api.MagnetismReflectometryReduction(InputWorkspace=self._event_workspace,
                                                 NormalizationWorkspace=ws_norm,
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
                                                 AngleOffset=angle_offset,
                                                 UseWLTimeAxis=False,
                                                 TimeAxisStep=self.configuration.tof_bins,
                                                 UseSANGLE=not self.configuration.use_dangle,
                                                 TimeAxisRange=self.tof_range,
                                                 SpecularPixel=self.configuration.peak_position,
                                                 ConstantQBinning=self.configuration.use_constant_q,
                                                 EntryName=str(self.entry_name),
                                                 OutputWorkspace=output_ws)

        ################## FOR COMPATIBILITY WITH QUICKNXS ##################
        run_object = ws.getRun()
        peak_min = run_object.getProperty("scatt_peak_min").value
        peak_max = run_object.getProperty("scatt_peak_max").value
        low_res_min = run_object.getProperty("scatt_low_res_min").value
        low_res_max = run_object.getProperty("scatt_low_res_max").value
        norm_x_min = run_object.getProperty("norm_peak_min").value
        norm_x_max = run_object.getProperty("norm_peak_max").value
        norm_y_min = run_object.getProperty("norm_low_res_min").value
        norm_y_max = run_object.getProperty("norm_low_res_max").value
        tth = ws.getRun().getProperty("SANGLE").getStatistics().mean * math.pi / 180.0
        quicknxs_scale = (float(norm_x_max)-float(norm_x_min)) * (float(norm_y_max)-float(norm_y_min))
        quicknxs_scale /= (float(peak_max)-float(peak_min)) * (float(low_res_max)-float(low_res_min))
        quicknxs_scale *= 0.005 / math.sin(tth)

        ws = api.Scale(InputWorkspace=output_ws, OutputWorkspace=output_ws,
                       factor=quicknxs_scale, Operation='Multiply')
        #####################################################################

        self.q = ws.readX(0)[:].copy()
        self._r = ws.readY(0)[:].copy() #* self.configuration.scaling_factor
        self._dr = ws.readE(0)[:].copy() #* self.configuration.scaling_factor

        #DeleteWorkspace(ws)
        self._reflectivity_workspace = str(ws)

    def offspec(self, direct_beam=None):
        """
            Extract off-specular scattering from 4D dataset (x,y,ToF,I).
            Uses a window in y to filter the 4D data
            and than sums all I values for each ToF and x channel.
            Qz,Qx,kiz,kfz is calculated using the x and ToF positions
            together with the tth-bank and direct pixel values.

            :param CrossSectionData direct_beam: if given, this data will be used to normalize the output
        """
        self.off_spec = off_specular.OffSpecular(self)
        return self.off_spec(direct_beam)

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

        phi = (np.arange(self.data.shape[1])[self.active_area_y[0]:
                                             self.active_area_y[1]]-y_pos)*rad_per_pixel

        v_edges = self.dist_mod_det/self.tof_edges*1e6 #m/s
        lambda_edges = H_OVER_M_NEUTRON/v_edges*1e10 #A
        wl = (lambda_edges[:-1] + lambda_edges[1:]) / 2.
        k = 2. * np.pi / wl

        # calculate ROI intensities and normalize by number of points
        P0 = self.configuration.cut_first_n_points
        PN = len(self.tof) - self.configuration.cut_last_n_points

        # calculate reciprocal space, incident and outgoing perpendicular wave vectors
        Qy = k[np.newaxis, np.newaxis, P0:PN]*(np.sin(phi)*np.cos(af)[:, np.newaxis])[:, :, np.newaxis]
        p_i = k[np.newaxis, np.newaxis, P0:PN]*((0*phi)+np.sin(ai)[:, np.newaxis])[:, :, np.newaxis]
        p_f = k[np.newaxis, np.newaxis, P0:PN]*((0*phi)+np.sin(af)[:, np.newaxis])[:, :, np.newaxis]
        Qz = p_i + p_f

        raw = self.data[self.active_area_x[0]:self.active_area_x[1],
                        self.active_area_y[0]:self.active_area_y[1],
                        P0:PN]

        intensity = scale * np.array(raw)
        d_intensity = scale * np.sqrt(raw)

        if direct_beam is not None:
            if not direct_beam.configuration.tof_bins == self.configuration.tof_bins:
                logging.error("Trying to normalize with a direct beam data set with different binning")

            norm_raw_multi_dim = direct_beam.data[self.active_area_x[0]:self.active_area_x[1],
                                                  self.active_area_y[0]:self.active_area_y[1], P0:PN]
            norm_raw = norm_raw_multi_dim.sum(axis=0).sum(axis=0)
            norm_d_raw = np.sqrt(norm_raw)

            surface = (self.active_area_x[1]-self.active_area_x[0]) * (self.active_area_y[1]-self.active_area_y[0])

            norm_raw /= surface * direct_beam.proton_charge * direct_beam.configuration.scaling_factor
            norm_d_raw /= surface * direct_beam.proton_charge * direct_beam.configuration.scaling_factor

            idxs = norm_raw > 0.
            d_intensity[:, :, idxs] = np.sqrt((d_intensity[:, :, idxs]/norm_raw[idxs][np.newaxis, np.newaxis, :])**2+
                                              (intensity[:, :, idxs]/norm_raw[idxs][np.newaxis, np.newaxis, :]**2*norm_d_raw[idxs][np.newaxis, np.newaxis, :])**2
                                             )
            intensity[:, :, idxs] /= norm_raw[idxs][np.newaxis, np.newaxis, :]
            intensity[:, :, np.logical_not(idxs)] = 0.
            d_intensity[:, :, np.logical_not(idxs)] = 0.

        # Create grid
        # bins=(self.options['gisans_gridy'], self.options['gisans_gridz']),
        #TODO: allow binning as application parameter
        self.SGrid, qy, qz = np.histogram2d(Qy.flatten(), Qz.flatten(),
                                            bins=(50, 50),
                                            weights=intensity.flatten())
        npoints, _, _ = np.histogram2d(Qy.flatten(), Qz.flatten(),
                                       bins=(50, 50))
        self.SGrid[npoints > 0] /= npoints[npoints > 0]
        self.SGrid = self.SGrid.transpose()
        qy = (qy[:-1]+qy[1:])/2.
        qz = (qz[:-1]+qz[1:])/2.
        self.QyGrid, self.QzGrid = np.meshgrid(qy, qz)

class NexusMetaData(object):
    """
        Class used to hold meta-data read before loading the neutron events
    """
    mid_q = 0
    is_direct_beam = False

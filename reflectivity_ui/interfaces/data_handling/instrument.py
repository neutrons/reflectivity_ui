"""
    This instrument description contains information
    that is instrument-specific and abstracts out how we obtain
    information from the data file
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long, bare-except
from __future__ import absolute_import, division, print_function
import sys
import math
import logging

# Import mantid according to the application configuration
from . import ApplicationConfiguration
application_conf = ApplicationConfiguration()
sys.path.insert(0, application_conf.mantid_path)
from mantid.simpleapi import *

# Option to use the slow flipper logs rather than the Analyzer/Polarizer logs
USE_SLOW_FLIPPER_LOG = True

# Constants
h = 6.626e-34  # m^2 kg s^-1
m = 1.675e-27  # kg


class Instrument(object):
    """
        Instrument class. Holds the data handling that is unique to a specific instrument.
    """
    n_x_pixel = 304
    n_y_pixel = 256
    huber_x_cut = 6.5
    peak_range_offset = 50
    tolerance = 0.05
    pixel_width = 0.0007
    instrument_name = "REF_M"
    instrument_dir = "/SNS/REF_M"
    file_search_template = "/SNS/REF_M/*/nexus/REF_M_%s"
    legacy_search_template = "/SNS/REF_M/*/data/REF_M_%s"

    # Filtering
    pol_state = 'SF1'
    pol_veto = 'SF1_Veto'
    ana_state = 'SF2'
    ana_veto = 'SF2_Veto'

    def __init__(self):
        logging.debug("Creating instrument")

    def dummy_filter_cross_sections(self, ws):
        """
            Filter events according to an aggregated state log.
            :param str file_path: file to read
    
            BL4A:SF:ICP:getDI
    
            015 (0000 1111): SF1=OFF, SF2=OFF, SF1Veto=OFF, SF2Veto=OFF
            047 (0010 1111): SF1=ON, SF2=OFF, SF1Veto=OFF, SF2Veto=OFF
            031 (0001 1111): SF1=OFF, SF2=ON, SF1Veto=OFF, SF2Veto=OFF
            063 (0011 1111): SF1=ON, SF2=ON, SF1Veto=OFF, SF2Veto=OFF
        """
        state_log = "BL4A:SF:ICP:getDI"
        states = {'Off_Off': 15,
                  'On_Off': 47,
                  'Off_On': 31,
                  'On_On': 63}
        cross_sections = []
    
        for pol_state in states:
            try:
                _ws = FilterByLogValue(InputWorkspace=ws, LogName=state_log, TimeTolerance=0.1,
                                      MinimumValue=states[pol_state],
                                      MaximumValue=states[pol_state], LogBoundary='Left',
                                      OutputWorkspace='%s_entry-%s' % (ws.getRunNumber(), pol_state))
                _ws.getRun()['cross_section_id'] = pol_state
                cross_sections.append(_ws)
            except:
                logging.error("Could not filter %s: %s", pol_state, sys.exc_info()[1])
    
        return cross_sections

    def load_data(self, file_path):
        """
            Load a data set according to the needs ot the instrument.
            Returns a WorkspaceGroup with any number of cross-sections.

            :param str file_path: path to the data file
        """
        if not USE_SLOW_FLIPPER_LOG:
            base_name = os.path.basename(file_path)
            _xs_list = MRFilterCrossSections(Filename=file_path,
                                            PolState=self.pol_state,
                                            AnaState=self.ana_state,
                                            PolVeto=self.pol_veto,
                                            AnaVeto=self.ana_veto,
                                            CrossSectionWorkspaces="%s_entry" % base_name)
            # Only keep good workspaced and get rid of the rejected events
            xs_list = [ws for ws in _xs_list if not ws.getRun()['cross_section_id'].value == 'unfiltered']
        else:
            ws = LoadEventNexus(Filename=file_path, OutputWorkspace="raw_events")
            xs_list = self.dummy_filter_cross_sections(ws)

        return xs_list

    @classmethod
    def scattering_angle(cls, ws, peak_position=None):
        """
            Determine the scattering angle in degrees
        """
        return MRGetTheta(ws, SpecularPixel=peak_position) * 180.0 / math.pi

    @classmethod
    def mid_q_value(cls, ws):
        """
            Get the mid q value, at the requested wl mid-point.
            :param workspace ws: Mantid workspace
        """
        wl = ws.getRun().getProperty('LambdaRequest').value[0]
        theta_d = Instrument.scattering_angle(ws) * math.pi / 180.0
        return 4.0*math.pi*math.sin(theta_d) / wl

    @classmethod
    def mid_q_value_from_data(cls, data_object):
        """
            Get the mid q value, at the requested wl mid-point.
            :param workspace ws: Mantid workspace
        """
        theta_d = (data_object.dangle - data_object.dangle0) / 2.0 * math.pi / 180.0
        return 4.0*math.pi*math.sin(theta_d) / data_object.lambda_center

    @classmethod
    def scattering_angle_from_data(cls, data_object):
        """
            Compute the scattering angle from a CrossSectionData object, in degrees.
            #TODO: this will go away once we keep the workspace

            @param data_object: CrossSectionData object
        """
        theta_d = (data_object.dangle - data_object.dangle0) / 2.0
        theta_d += ((data_object.dpix - data_object.configuration.peak_position) * cls.pixel_width) * 180.0 / math.pi / (2.0 * data_object.dist_sam_det)
        return theta_d

    def check_direct_beam(self, ws, peak_position=None):
        """
            Determine whether this data is a direct beam
        """
        sangle = ws.getRun().getProperty("SANGLE").getStatistics().mean
        theta = self.scattering_angle(ws, peak_position)
        huber_x = ws.getRun().getProperty("HuberX").getStatistics().mean
        return not ((theta > self.tolerance or sangle > self.tolerance) and huber_x < self.huber_x_cut)

    def direct_beam_match(self, scattering, direct_beam, skip_slits=False):
        """
            Verify whether two data sets are compatible.
        """
        if scattering.number == direct_beam.number \
            or ((direct_beam.scattering_angle > self.tolerance \
                 or direct_beam.sangle > self.tolerance) and direct_beam.huber_x < self.huber_x_cut):
            logging.error("Run %s may not be a direct beam", direct_beam.number)

        if math.fabs(scattering.lambda_center-direct_beam.lambda_center) < self.tolerance \
            and (skip_slits or \
            (math.fabs(scattering.slit1_width-direct_beam.slit1_width) < self.tolerance \
            and math.fabs(scattering.slit2_width-direct_beam.slit2_width) < self.tolerance \
            and math.fabs(scattering.slit3_width-direct_beam.slit3_width) < self.tolerance)):
            return True
        return False

    @classmethod
    def get_info(cls, workspace, data_object):
        """
            Retrieve information that is specific to this particular instrument

            @param workspace: Mantid workspace
            @param data_object: CrossSectionData object
        """
        data = workspace.getRun()
        data_object.lambda_center=data['LambdaRequest'].value[0]
        data_object.dangle=data['DANGLE'].value[0]
        data_object.dangle0=data['DANGLE0'].value[0]
        data_object.dpix=data['DIRPIX'].value[0]
        data_object.slit1_width=data['S1HWidth'].value[0]
        data_object.slit2_width=data['S2HWidth'].value[0]
        data_object.slit3_width=data['S3HWidth'].value[0]
        data_object.huber_x=data['HuberX'].getStatistics().mean

        data_object.sangle=data['SANGLE'].value[0]

        data_object.dist_sam_det=data['SampleDetDis'].value[0]*1e-3
        data_object.dist_mod_det=data['ModeratorSamDis'].value[0]*1e-3+data_object.dist_sam_det
        data_object.dist_mod_mon=data['ModeratorSamDis'].value[0]*1e-3-2.75

        # Get these from instrument
        data_object.pixel_width = float(workspace.getInstrument().getNumberParameter("pixel-width")[0]) / 1000.0
        data_object.n_det_size_x = int(workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0]) # 304
        data_object.n_det_size_y = int(workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0]) # 256
        data_object.det_size_x = data_object.n_det_size_x * data_object.pixel_width # horizontal size of detector [m]
        data_object.det_size_y = data_object.n_det_size_y * data_object.pixel_width # vertical size of detector [m]

        # The following active area used to be taken from instrument.DETECTOR_REGION
        data_object.active_area_x = (8, 295)
        data_object.active_area_y = (8, 246)

        # Convert to standard names
        data_object.direct_pixel = data_object.dpix
        data_object.angle_offset = data_object.dangle0

    def integrate_detector(self, ws, specular=True):
        """
            Integrate a workspace along either the main direction (specular=False) or
            the low-resolution direction (specular=True.

            :param ws: Mantid workspace
            :param specular bool: if True, the low-resolution direction is integrated over
        """
        ws_summed = RefRoi(InputWorkspace=ws, IntegrateY=specular,
                           NXPixel=self.n_x_pixel, NYPixel=self.n_y_pixel,
                           ConvertToQ=False,
                           OutputWorkspace="ws_summed")

        integrated = Integration(ws_summed)
        integrated = Transpose(integrated)
        return integrated

    @classmethod
    def find_direct_beam(cls):
        """
            Find a direct beam suitable for this data set by looking into
            its data directory.
            :param ws: Mantid workspace
        """
        
        return
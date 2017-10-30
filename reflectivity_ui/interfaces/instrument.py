"""
    This instrument description contains information
    that is instrument-specific and abstracts out how we obtain
    information from the data file
"""
#pylint: disable=invalid-name, too-many-instance-attributes, line-too-long
from __future__ import absolute_import, division, print_function

class Instrument(object):
    """
        Instrument class. Holds the data handling that is unique to a specific instrument.
    """
    def __init__(self):
        self.instrument_name = "REFM"
        self.tof_range = [0,0]

    def get_tof_range(self, run_object):
        """
            Determine TOF range from the data
        """
        sample_detector_distance = run_object['SampleDetDis'].getStatistics().mean / 1000.0
        source_sample_distance = run_object['ModeratorSamDis'].getStatistics().mean / 1000.0
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

        #TODO: these don't exist in the DASLogs
        #data_object.slit1_dist=-data['instrument/aperture1/distance'].value[0]*1000.
        #data_object.slit2_dist=-data['instrument/aperture2/distance'].value[0]*1000.
        #data_object.slit3_dist=-data['instrument/aperture3/distance'].value[0]*1000.

        data_object.sangle=data['SANGLE'].value[0]

        data_object.dist_sam_det=data['SampleDetDis'].value[0]*1e-3
        data_object.dist_mod_det=data['ModeratorSamDis'].value[0]*1e-3+data_object.dist_sam_det
        data_object.dist_mod_mon=data['ModeratorSamDis'].value[0]*1e-3-2.75

        # Get these from instrument
        data_object.det_size_x = int(workspace.getInstrument().getNumberParameter("number-of-x-pixels")[0]) #304
        data_object.det_size_y = int(workspace.getInstrument().getNumberParameter("number-of-y-pixels")[0]) #256

        # The following active area used to be taken from instrument.DETECTOR_REGION
        data_object.active_area_x = (8, 295)
        data_object.active_area_y = (8, 246)

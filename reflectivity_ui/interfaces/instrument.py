"""
    This instrument description contains information
    that is instrument-specific and abstracts out how we obtain
    information from the data file
"""


class Instrument(object):
    def __init__(self):
        self.instrument_name = "REFM"

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
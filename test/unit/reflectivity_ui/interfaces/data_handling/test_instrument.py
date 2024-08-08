# package imports
from reflectivity_ui.interfaces.configuration import Configuration

# 3rd party imports
import pytest


@pytest.mark.datarepo
def test_load_data_deadtime(data_server):
    """Test load data with and without dead-time correction"""
    conf = Configuration()
    file_path = data_server.path_to("REF_M_42112")
    corrected_events = [52283.51, 42028.15, 66880.96, 43405.89]

    # load with dead-time correction
    conf.apply_deadtime = True
    ws_list = conf.instrument.load_data(file_path, conf)
    assert len(ws_list) == 4
    for iws, ws in enumerate(ws_list):
        assert "dead_time_applied" in ws.getRun()
        assert ws.extractY().sum() == pytest.approx(corrected_events[iws])

    # load without dead-time correction
    conf.apply_deadtime = False
    ws_list = conf.instrument.load_data(file_path, conf)
    assert len(ws_list) == 4
    for ws in ws_list:
        assert "dead_time_applied" not in ws.getRun()
        assert ws.extractY().sum() == ws.getNumberEvents()

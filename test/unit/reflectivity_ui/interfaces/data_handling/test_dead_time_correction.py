# package imports
from reflectivity_ui.interfaces.data_handling.DeadTimeCorrection import SingleReadoutDeadTimeCorrection

# 3rd-party imports
import mantid.simpleapi as api
import pytest
from mantid.kernel import amend_config

# standard imports


@pytest.mark.datarepo
@pytest.mark.parametrize("is_paralyzable", [False, True])
def test_deadtime(is_paralyzable, data_server):
    """Test of the dead-time correction algorithm SingleReadoutDeadTimeCorrection"""
    with amend_config(data_dir=data_server.h5_full_path):
        ws = api.Load("REF_M_42112")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.setProperty("Paralyzable", is_paralyzable)
    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)
    for c in corr:
        assert 1.0 <= c < 1.001, "value not between 1.0 and 1.001"

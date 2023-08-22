# package imports
from reflectivity_ui.interfaces.configuration import Configuration
from reflectivity_ui.interfaces.data_handling.data_manipulation import smart_stitch_reflectivity
from reflectivity_ui.interfaces.data_handling.data_set import NexusData
from reflectivity_ui.interfaces.data_manager import DataManager

# 3rd-party imports
import copy
import mantid.simpleapi as api
import numpy as np
import pytest
from test import SNS_REFM_MOUNTED


mock_reduced_file_str = (
    "# Datafile created by QuickNXS 3.1.0.dev2\n"
    "# Datafile created using Mantid 6.6.0\n"
    "# Date: 2023-07-05 13:31:21\n"
    "# Type: Specular\n"
    "# Input file indices: 42112,42116,42113\n"
    "# Extracted states: -+\n"
    "#\n"
    "# [Direct Beam Runs]\n"
    "#    DB_ID        P0        PN     x_pos   x_width     y_pos   y_width    bg_pos  bg_width      dpix       tth    number      File\n"  # noqa: E501
    "#        1         0         0       195        12     126.5       155        55        24       194         0     42099  /SNS/REF_M/IPTS-30794/nexus/REF_M_42099.nxs.h5\n"  # noqa: E501
    "#        2         0         0       195        12       127       154        55        24       194         0     42100  /SNS/REF_M/IPTS-30794/nexus/REF_M_42100.nxs.h5\n"  # noqa: E501
    "#        3         0         0       195        12       127       154        55        24       194         0     42100  /SNS/REF_M/IPTS-30794/nexus/REF_M_42100.nxs.h5\n"  # noqa: E501
    "#\n"
    "# [Data Runs]\n"
    "#    scale        P0        PN     x_pos   x_width     y_pos   y_width    bg_pos  bg_width       fan      dpix       tth    number     DB_ID      File\n"  # noqa: E501
    "#        1        15        10       167        12     163.5      72.9        55        24     False       194  0.00653668     42112         1  /SNS/REF_M/IPTS-30794/nexus/REF_M_42112.nxs.h5\n"  # noqa: E501
    "# 0.183654        15        10     189.3        12     162.2      68.3        55        24     False       194  0.799577     42116         2  /SNS/REF_M/IPTS-30794/nexus/REF_M_42116.nxs.h5\n"  # noqa: E501
    "#   0.1375        15        10     167.5        12     159.9      67.4        55        24     False       194  0.798876     42113         3  /SNS/REF_M/IPTS-30794/nexus/REF_M_42113.nxs.h5\n"  # noqa: E501
    "#\n"
    "# [Global Options]\n"
    "# name           value\n"
    "# sample_length  10.0\n"
    "#\n"
    "# [Data]\n"
    "#     Qz [1/A]	    R [a.u.]	   dR [a.u.]	   dQz [1/A]"
)


@pytest.fixture
def mocker_file_open(mocker):
    mocked_reduced_file = mocker.mock_open(read_data=mock_reduced_file_str)
    mocker.patch("builtins.open", mocked_reduced_file)


@pytest.fixture
def stitching_config():
    settings = {
        "use_dangle": True,
        "normalize_to_unity": False,
        "total_reflectivity_q_cutoff": 0.008,
        "do_final_rebin": False,
        "final_rebin_step": -0.01,
        "match_direct_beam": True,
    }
    config = Configuration()
    for key in settings:
        setattr(config, key, settings[key])
    return config


@pytest.fixture
def stitching_reduction_list():
    """List of NexusData objects for testing of stitching"""

    class _MockCrossSectionData(object):
        """Test class to use instead of CrossSectionData, which requires EventWorkspaces"""

        def __init__(self, xs: str, config: Configuration, data_x: list, data_y: list, ws_name: str):
            self.name = xs
            self.configuration = copy.deepcopy(config)
            data_e = [0.1 * y for y in data_y]
            workspace = api.CreateWorkspace(DataX=data_x, DataY=data_y, DataE=data_e, OutputWorkspace=ws_name)
            self.ws = workspace
            self.q = self.ws.readX(0)[:].copy()
            self._r = np.ma.masked_equal(self.ws.readY(0)[:].copy(), 0)
            self._dr = self.ws.readE(0)[:].copy()

    configuration = Configuration()
    configuration.cut_first_n_points = 0
    configuration.cut_last_n_points = 0

    # create curve 1
    data_x = [1.0, 2.0, 3.0, 4.0]
    xs1_on_on = _MockCrossSectionData("On_On", configuration, data_x, [5.0, 5.0, 5.0, 5.0], "run1_on_on")
    xs1_on_off = _MockCrossSectionData("On_Off", configuration, data_x, [7.0, 7.0, 7.0, 7.0], "run1_on_off")
    curve1 = NexusData("path", configuration)
    curve1.cross_sections = {"On_On": xs1_on_on, "On_Off": xs1_on_off}

    # create curve 2
    data_x = [4.0, 5.0, 6.0, 7.0]
    xs2_on_on = _MockCrossSectionData("On_On", configuration, data_x, [3.0, 3.0, 3.0, 3.0], "run2_on_on")
    xs2_on_off = _MockCrossSectionData("On_Off", configuration, data_x, [2.0, 2.0, 2.0, 2.0], "run2_on_off")
    curve2 = NexusData("path", configuration)
    curve2.cross_sections = {"On_On": xs2_on_on, "On_Off": xs2_on_off}

    # create curve 3
    data_x = [7.0, 8.0, 9.0, 10.0]
    xs3_on_on = _MockCrossSectionData("On_On", configuration, data_x, [1.0, 1.0, 1.0, 1.0], "run3_on_on")
    xs3_on_off = _MockCrossSectionData("On_Off", configuration, data_x, [1.0, 1.0, 1.0, 1.0], "run3_on_off")
    curve3 = NexusData("path", configuration)
    curve3.cross_sections = {"On_On": xs3_on_on, "On_Off": xs3_on_off}

    return [curve1, curve2, curve3]


class TestDataManipulation(object):
    @pytest.mark.skipif(not SNS_REFM_MOUNTED, reason="/SNS/REF_M/ is not mounted")
    def test_smart_stitch_reflectivity(self, data_server, mocker_file_open, stitching_config):
        manager = DataManager(data_server.directory)
        manager.load_data_from_reduced_file("note: file read is mocked", stitching_config)
        if len(manager.reduction_list) < 1:
            raise IOError("Files missing.")
        scaling_factors, scaling_errors = smart_stitch_reflectivity(manager.reduction_list, "Off_On", False, 0.008)
        assert scaling_factors == pytest.approx([1.0, 0.1809, 0.1556], abs=0.001)
        assert scaling_errors == pytest.approx([0.0, 0.0, 0.0], abs=0.001)

    @pytest.mark.parametrize(
        "normalize_to_unity, q_cutoff, global_fit, expected_scaling_factors, expected_scaling_errors",
        [
            (False, None, False, [1.0, 1.66667, 5.0], [0.0, 0.23570, 1.00000]),
            (True, 1.5, False, [0.2, 0.33333, 1.0], [0.02, 0.057735, 0.22361]),
            (False, None, True, [1.0, 2.4, 6.0], [0.0, 0.24403, 0.85988]),
            (True, 1.5, True, [0.2, 0.48, 1.2], [0.02, 0.068455, 0.20970]),
        ],
    )
    def test_smart_stitch_parameters(
        self,
        stitching_reduction_list,
        normalize_to_unity,
        q_cutoff,
        global_fit,
        expected_scaling_factors,
        expected_scaling_errors,
    ):
        """Test all combinations of the smart_stitch_reflectivity parameters `normalize_to_unity` and `global_fit`

        The fixture stitching_reduction_list has three runs with two cross-sections each: `On_On` and `On_Off`.
        When `global_fit` is True, both cross-sections are used to calculate the scaling factors.
        """
        scaling_factors, scaling_errors = smart_stitch_reflectivity(
            stitching_reduction_list, "On_On", normalize_to_unity, q_cutoff, global_fit
        )
        assert scaling_factors == pytest.approx(expected_scaling_factors, abs=0.001)
        assert scaling_errors == pytest.approx(expected_scaling_errors, abs=0.001)
        # Delete workspaces
        api.mtd.clear()


if __name__ == "__main__":
    pytest.main([__file__])

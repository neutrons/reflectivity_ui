import pytest

from reflectivity_ui.interfaces.configuration import Configuration
from reflectivity_ui.interfaces.data_handling.data_manipulation import smart_stitch_reflectivity
from reflectivity_ui.interfaces.data_manager import DataManager

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


class TestDataManipulation(object):

    @pytest.mark.skipif(not SNS_REFM_MOUNTED, reason="/SNS/REF_M/ is not mounted")
    def test_smart_stitch_reflectivity(self, data_server, mocker_file_open, stitching_config):
        manager = DataManager(data_server.directory)
        manager.load_data_from_reduced_file("note: file read is mocked", stitching_config)
        if len(manager.reduction_list) < 1:
            raise IOError("Files missing.")
        scaling_factors = smart_stitch_reflectivity(manager.reduction_list, None, False, 0.008)
        assert scaling_factors == pytest.approx([1.0, 0.1809, 0.1354], abs=0.001)


if __name__ == "__main__":
    pytest.main([__file__])

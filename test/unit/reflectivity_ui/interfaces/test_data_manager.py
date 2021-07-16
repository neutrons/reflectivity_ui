from __future__ import absolute_import, division, print_function, unicode_literals

# local imports
from reflectivity_ui.interfaces.data_manager import DataManager
from reflectivity_ui.interfaces.configuration import Configuration
import reflectivity_ui.interfaces.data_handling.data_manipulation as dm

# 3rd-party imports
import pytest


class TestDataManagerTest(object):

    def test_manager(self, data_server):
        manager = DataManager(data_server.directory)
        manager.load(data_server.file("REF_M_29160"), Configuration())

        assert manager.current_file == "REF_M_29160"

        manager.add_active_to_reduction()
        assert manager.find_data_in_reduction_list(manager._nexus_data) == 0
        assert manager.find_data_in_direct_beam_list(manager._nexus_data) is None

        q_range = manager._nexus_data.get_q_range()
        assert q_range[0: 2] == pytest.approx([0.034, 0.068], abs=0.05)
        assert manager.add_active_to_normalization()
        assert manager.remove_active_from_normalization() == 0

        manager.set_active_data_from_reduction_list(0)
        manager.set_active_data_from_direct_beam_list(0)
        manager.calculate_reflectivity()
        manager.calculate_reflectivity(specular=False)
        manager.strip_overlap()

        dm.generate_script(manager.reduction_list, manager.reduction_states[0])
        dm.stitch_reflectivity(manager.reduction_list)
        dm.merge_reflectivity(manager.reduction_list, manager.reduction_states[0])
        dm.get_scaled_workspaces(manager.reduction_list, manager.reduction_states[0])
        dm.stitch_reflectivity(manager.reduction_list)

    def test_load_reduced(self, data_server):
        manager = DataManager(data_server.directory)
        manager.load_data_from_reduced_file(data_server.file('REF_M_29160_Specular_++.dat'))


if __name__ == '__main__':
    pytest.main([__file__])

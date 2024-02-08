# local imports
from reflectivity_ui.interfaces.configuration import Configuration
from reflectivity_ui.interfaces.data_handling.data_set import CrossSectionData, NexusData
from reflectivity_ui.interfaces.main_window import MainWindow
from test import SNS_REFM_MOUNTED
from test.ui import ui_utilities

# third party imports
import pytest
from qtpy import QtCore, QtWidgets


# standard library imports


class TestMainGui:
    def test_init(self, qtbot):
        window_main = MainWindow()
        qtbot.addWidget(window_main)
        assert "QuickNXS Magnetic Reflectivity" in window_main.windowTitle()

    @pytest.mark.skipif(not SNS_REFM_MOUNTED, reason="/SNS/REF_M/ is not mounted")
    def test_enter_run_number(self, qtbot):
        window_main = MainWindow()
        qtbot.addWidget(window_main)
        ui_utilities.setText(window_main.numberSearchEntry, "42100", press_enter=True)

    def test_set_global_stitching(self, qtbot):
        """Test that the configuration is updated based on the checkbox value"""
        window_main = MainWindow()
        qtbot.addWidget(window_main)
        window_main.global_fit_checkbox.setChecked(True)
        assert window_main.file_handler.get_configuration().global_stitching is True
        window_main.global_fit_checkbox.setChecked(False)
        assert window_main.file_handler.get_configuration().global_stitching is False

    def test_active_channel(self, mocker, qtbot):
        """Test that selecting a cross-section radio button updates the active channel"""
        # mock updating the plots
        mocker.patch("reflectivity_ui.interfaces.main_window.MainWindow.plotActiveTab", return_value=True)
        mocker.patch("reflectivity_ui.interfaces.plotting.PlotManager.plot_refl", return_value=True)
        mocker.patch("reflectivity_ui.interfaces.plotting.PlotManager.plot_projections", return_value=True)

        # create the SUT
        window_main = MainWindow()
        qtbot.addWidget(window_main)

        # set up data objects for two channels
        configuration = window_main.file_handler.get_configuration()
        channel0 = CrossSectionData("On_Off", configuration)
        channel1 = CrossSectionData("On_On", configuration)
        nexus_data = NexusData("filepath", configuration)
        nexus_data.cross_sections = {channel0.name: channel0, channel1.name: channel1}
        window_main.data_manager._nexus_data = nexus_data
        window_main.data_manager.set_channel(0)

        assert window_main.data_manager.active_channel.name == channel0.name

        # change the selected channel
        window_main.selectedChannel1.setChecked(True)

        # check the active channel in the data manager
        assert window_main.data_manager.active_channel.name == channel1.name
        # check the current channel name displayed in the UI
        assert channel1.name in window_main.ui.currentChannel.text()

    @pytest.mark.parametrize("table_widget", ["reductionTable", "normalizeTable"])
    def test_reduction_table_right_click(self, table_widget, qtbot, mocker):
        mock_save_run_data = mocker.patch(
            "reflectivity_ui.interfaces.event_handlers.main_handler.MainHandler.save_run_data"
        )
        window_main = MainWindow()
        qtbot.addWidget(window_main)
        window_main.data_manager.reduction_list = [NexusData("filepath", Configuration())]
        window_main.data_manager.direct_beam_list = [NexusData("filepath", Configuration())]
        table = getattr(window_main.ui, table_widget)
        table.insertRow(0)

        def handle_menu():
            """Press Enter on item in menu and check that the function was called"""
            menu = table.findChild(QtWidgets.QMenu)
            action = menu.actions()[0]
            assert action.text() == "Export data"
            qtbot.keyClick(menu, QtCore.Qt.Key_Down)
            qtbot.keyClick(menu, QtCore.Qt.Key_Enter)
            mock_save_run_data.assert_called_once()

        QtCore.QTimer.singleShot(200, handle_menu)
        pos = QtCore.QPoint()
        table.customContextMenuRequested.emit(pos)


if __name__ == "__main__":
    pytest.main([__file__])

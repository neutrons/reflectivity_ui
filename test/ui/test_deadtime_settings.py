# local imports
from reflectivity_ui.interfaces.configuration import Configuration
from reflectivity_ui.interfaces.main_window import MainWindow

# standard imports
import unittest.mock as mock

# third-party imports
import pytest
from qtpy import QtCore, QtWidgets


def test_show_deadtime_settings_default_values(qtbot):
    """Test showing the deadtime settings dialog and keeping the default values"""
    main_window = MainWindow()
    qtbot.addWidget(main_window)

    def display_deadtime_settings():
        # open the deadtime settings dialog and close without changing values
        dialog = QtWidgets.QApplication.activeModalWidget()
        # press Enter to accept (click Ok)
        qtbot.keyClick(dialog, QtCore.Qt.Key_Enter, delay=1)

    QtCore.QTimer.singleShot(200, display_deadtime_settings)
    main_window.open_deadtime_settings()

    assert Configuration.paralyzable_deadtime == True
    assert Configuration.deadtime_value == 4.2
    assert Configuration.deadtime_tof_step == 100


def test_show_deadtime_settings_updated_values(qtbot):
    """Test showing the deadtime settings dialog and updating the values"""
    new_paralyzable = False
    new_dead_time = 5.0
    new_tof_step = 200

    main_window = MainWindow()
    qtbot.addWidget(main_window)

    def update_deadtime_settings():
        # update the values in the deadtime settings dialog and close the dialog
        dialog = QtWidgets.QApplication.activeModalWidget()
        paralyzable_checkbox = dialog.findChild(QtWidgets.QCheckBox)
        paralyzable_checkbox.setChecked(new_paralyzable)
        deadtime_spinbox = dialog.findChild(QtWidgets.QDoubleSpinBox, "dead_time_value")
        deadtime_spinbox.setValue(new_dead_time)
        tof_spinbox = dialog.findChild(QtWidgets.QDoubleSpinBox, "dead_time_tof")
        tof_spinbox.setValue(new_tof_step)
        # press Enter to accept (click Ok)
        qtbot.keyClick(dialog, QtCore.Qt.Key_Enter, delay=1)

    QtCore.QTimer.singleShot(200, update_deadtime_settings)
    main_window.open_deadtime_settings()

    assert Configuration.paralyzable_deadtime == new_paralyzable
    assert Configuration.deadtime_value == new_dead_time
    assert Configuration.deadtime_tof_step == new_tof_step


if __name__ == "__main__":
    pytest.main([__file__])

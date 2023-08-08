# local imports
from reflectivity_ui.interfaces.main_window import MainWindow
from reflectivity_ui.interfaces.data_handling.data_set import REFLECTIVITY_THRESHOLD_VALUE
from test import SNS_REFM_MOUNTED
from test.ui import ui_utilities

# third party imports
import numpy as np
from PyQt5 import QtCore
import pytest

# standard library imports


@pytest.mark.skipif(not SNS_REFM_MOUNTED, reason="/SNS/REF_M/ is not mounted")
def test_missing_cross_section(qtbot):
    r"""Test a run where the cross section corresponding to the On-On spin combination has no integrated
    proton charge. The application produces and empty reflectivity curve for On-On."""
    main_window = MainWindow()
    qtbot.addWidget(main_window)
    # load the run and find the total "intensity" of the x vs TOF plot
    ui_utilities.setText(main_window.numberSearchEntry, "42100", press_enter=True)
    intensity_off_on = np.sum(ui_utilities.data_from_plot2D(main_window.xtof_overview))
    # select the On-On spin combination
    qtbot.mouseClick(main_window.selectedChannel1, QtCore.Qt.LeftButton)
    # check that the reflectivity curve is empty
    _, data_y = ui_utilities.data_from_plot1D(main_window.refl)
    assert np.all(data_y <= REFLECTIVITY_THRESHOLD_VALUE)
    # check the x vs TOF plot has changed
    intensity_on_on = np.sum(ui_utilities.data_from_plot2D(main_window.xtof_overview))
    np.testing.assert_almost_equal(intensity_off_on / intensity_on_on, 279, decimal=0)


if __name__ == '__main__':
    pytest.main([__file__])

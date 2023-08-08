# local imports
from reflectivity_ui.interfaces.main_window import MainWindow
from test import SNS_REFM_MOUNTED
from test.ui import ui_utilities

# third party imports
import pytest

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


if __name__ == "__main__":
    pytest.main([__file__])

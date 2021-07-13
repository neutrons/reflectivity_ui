from reflectivity_ui.interfaces.event_handlers.widgets import CustomDialog
from PyQt5 import QtCore

wait = 100


def test_customized_dialog(qtbot):
    """Test customized dialog box
    """
    window = CustomDialog(None, 'pyqt test', 'Hello World')

    qtbot.addWidget(window)
    window.show()
    qtbot.wait(wait)

    qtbot.mouseClick(window.buttonBox, QtCore.Qt.LeftButton)
    qtbot.wait(wait)

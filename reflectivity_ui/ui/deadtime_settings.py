# third-party imports
from qtpy.QtWidgets import QDialog, QWidget
from xml.dom.minidom import Document, Element
from typing import Any, Callable, Dict
import os

from qtpy.uic import loadUi


class DeadTimeSettingsModel():
    """Stores all options for the dead time correction. These are global options"""

    apply_deadtime: bool = False
    paralyzable: bool = True
    dead_time: float = 4.2
    tof_step: float = 100.0

class DeadTimeSettingsView(QDialog):
    """
    Dialog to choose the dead time correction options.
    """

    def __init__(self, parent: QWidget):
        super().__init__(parent)
        filepath = os.path.join(os.path.dirname(__file__), "deadtime_settings.ui")
        self.ui = loadUi(filepath, baseinstance=self)
        self.options = self.get_state_from_form()

    def set_state(self, paralyzable, dead_time, tof_step):
        """
        Store options and populate the form
        :param apply_correction: If True, dead time correction will be applied
        :param paralyzable: If True, a paralyzable correction will be used
        :param dead_time: Value of the dead time in micro second
        :param tof_step: TOF binning in micro second
        """
        self.ui.use_paralyzable.setChecked(paralyzable)
        self.ui.dead_time_value.setValue(dead_time)
        self.ui.dead_time_tof.setValue(tof_step)
        self.options = self.get_state_from_form()

    def get_state_from_form(self) -> dict:
        r"""Read the options from the form.

        Returns
        -------
        Dictionary whose keys must match fields of class `DeadTimeSettingsModel`
        """
        return {
            'paralyzable': self.ui.use_paralyzable.isChecked(),
            'dead_time': self.ui.dead_time_value.value(),
            'tof_step': self.ui.dead_time_tof.value(),
        }

    def accept(self):
        """
        Read in the options on the form when the OK button is
        clicked.
        """
        self.options = self.get_state_from_form()
        self.close()

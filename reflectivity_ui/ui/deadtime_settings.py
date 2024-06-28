# package imports
from reflectivity_ui.interfaces.configuration import Configuration

# third-party imports
from qtpy.QtWidgets import QDialog, QWidget
from qtpy.uic import loadUi

# standard imports
import os


class DeadTimeSettingsView(QDialog):
    """
    Dialog to choose the dead time correction options.
    """

    def __init__(self, parent: QWidget):
        super().__init__(parent)
        filepath = os.path.join(os.path.dirname(__file__), "deadtime_settings.ui")
        self.ui = loadUi(filepath, baseinstance=self)
        self.set_state_from_global_config()

    def set_state_from_global_config(self):
        """
        Populate the form with the current global configuration
        """
        self.ui.use_paralyzable.setChecked(Configuration.paralyzable_deadtime)
        self.ui.dead_time_value.setValue(Configuration.deadtime_value)
        self.ui.dead_time_tof.setValue(Configuration.deadtime_tof_step)

    def accept(self):
        """
        Read in the options on the form when the OK button is
        clicked and update the global configuration.
        """
        Configuration.paralyzable_deadtime = self.ui.use_paralyzable.isChecked()
        Configuration.deadtime_value = self.ui.dead_time_value.value()
        Configuration.deadtime_tof_step = self.ui.dead_time_tof.value()
        self.close()

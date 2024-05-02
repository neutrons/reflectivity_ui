from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QCheckBox, QDoubleSpinBox, QSpinBox

from reflectivity_ui.interfaces.configuration import Configuration


class ConfigurationHandler:
    """
    Handles updating the configuration state upon changes in the UI
    """

    def __init__(self, main_window):
        self.main_window = main_window
        self.ui = main_window.ui
        self.connect_config_events()

    def config_setter_factory(self, config_name: str, is_checkbox: bool):
        """
        Generate anonymous functions to serve as callback when any of the global configurations
        (`Configuration` class variables) are updated in the UI.

        Each callback will be associated to one configuration parameter. Upon invoked,
        the `Configuration` class variable value will be updated.

        Parameters
        ----------
        config_name: str
            Name of the Configuration variable to update
        is_checkbox: bool
            True if the widget type is QCheckBox
        """

        def config_setter(value: int | float):
            if is_checkbox:
                bool_value = value == Qt.CheckState.Checked
                setattr(Configuration, config_name, bool_value)
            else:
                setattr(Configuration, config_name, value)

        return config_setter

    def connect_config_events(self):
        widget_names = [
            "final_rebin_checkbox",
            "q_rebin_spinbox",
            "normalize_to_unity_checkbox",
            "normalization_q_cutoff_spinbox",
            "global_fit_checkbox",
            "polynomial_stitching_degree_spinbox",
            "polynomial_stitching_points_spinbox",
            "polynomial_stitching_checkbox",
            "fanReflectivity",
            "sample_size_spinbox",
            "bandwidth_spinbox",
        ]
        config_names = [
            "do_final_rebin",
            "final_rebin_step",
            "normalize_to_unity",
            "total_reflectivity_q_cutoff",
            "global_stitching",
            "polynomial_stitching_degree",
            "polynomial_stitching_points",
            "polynomial_stitching",
            "use_constant_q",
            "sample_size",
            "wl_bandwidth",
        ]

        for widget_name, config_name in zip(widget_names, config_names):
            # get the widget signal to connect
            widget = getattr(self.ui, widget_name)
            is_checkbox = False
            if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
                signal_name = "valueChanged"
            elif isinstance(widget, QCheckBox):
                is_checkbox = True
                signal_name = "stateChanged"
            else:
                raise ValueError(f"{type(widget)} not supported by ConfigurationHandler")
            signal = getattr(widget, signal_name)
            # connect the signal to the updater
            signal.connect(self.config_setter_factory(config_name, is_checkbox))

# local imports


# third party imports
from numpy.ma import MaskedArray
from PyQt5 import QtCore

# standard library imports


def setText(widget, text, press_enter=True):
    r"""Set the text of a widget and optionally press enter."""
    assert getattr(widget, "setText", None) is not None
    widget.setText(text)
    if press_enter:
        widget.setFocus(QtCore.Qt.MouseFocusReason)
        widget.returnPressed.emit()


def data_from_plot1D(widget: "MplWidget", line_number=0) -> tuple:
    r"""Get the data from an MplWidget representing a 1D plot
    Returns
    -------
    X and Y data as a tuple of numpy arrays
    """
    figure = widget.canvas.fig
    axes = figure.get_axes()[0]
    return axes.get_lines()[line_number].get_data()


def data_from_plot2D(widget: "MplWidget") -> MaskedArray:
    r"""Get the data from an MplWidget representing a 1D plot
    Returns
    -------
    2D data as a masked numpy array
    """
    figure = widget.canvas.fig
    axes = figure.get_axes()[0]
    return axes.get_images()[0].get_array()

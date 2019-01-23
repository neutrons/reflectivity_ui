"""
   Dialog to show final reduced data.
"""
#pylint: disable=bare-except
from __future__ import absolute_import, division, print_function, unicode_literals
from PyQt5 import QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_result_viewer

class ResultViewer(QtWidgets.QDialog, reflectivity_ui.interfaces.generated.ui_result_viewer.Ui_Dialog):
    """
        Reduction dialog
    """
    default_template = u'(instrument)_{numbers}_{item}_{state}.{type}'

    def __init__(self, parent):
        super(ResultViewer, self).__init__(parent)
        self.setupUi(self)
        self.resize(1024, 1024)
        self.settings = QtCore.QSettings('.refredm')

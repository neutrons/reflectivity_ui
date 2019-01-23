"""
   Dialog to show final reduced data.
"""
#pylint: disable=bare-except
from __future__ import absolute_import, division, print_function, unicode_literals
import os
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

        self.settings = QtCore.QSettings('.refredm')


    def accept(self):
        """
            Save the current options and close dialog
        """
        #self.save_settings()
        #self.is_accepted = True
        self.close()

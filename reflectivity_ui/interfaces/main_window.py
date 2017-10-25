#pylint: disable=invalid-name, line-too-long
"""
    Main application window
"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import logging
import glob
from PyQt4 import QtGui, QtCore
import reflectivity_ui.interfaces.generated.ui_main_window

class MainWindow(QtGui.QMainWindow, reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow):
    """
        Main applicqtion window
    """
    def __init__(self):
        """
            Initialization
        """
        # Base class
        QtGui.QMainWindow.__init__(self)

        # Initialize the UI widgets
        self.ui = reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow()
        self.ui.setupUi(self)

        # Application settings
        self.settings = QtCore.QSettings('.refredm')
        self._current_directory = self.settings.value('current_directory', os.path.expanduser('~')).toString()
        self._current_file = None

        # Initialization for specific instrument
        # Retrieve configuration from config and enable/disable features

        self.hide_unsupported()

    def hide_unsupported(self):
        """
            Hide what we don't support
        """
        # Hide event filtering (which is not really event filtering)
        for i in range(self.ui.event_filtering_layout.rowCount()*self.ui.event_filtering_layout.columnCount()):
            if self.ui.event_filtering_layout.itemAt(i):
                self.ui.event_filtering_layout.itemAt(i).widget().hide()

        # Hide format selectiuon since we're only using event data
        self.ui.oldFormatActive.hide()
        self.ui.histogramActive.hide()
        self.ui.eventActive.hide()

        # Hide quick reduce button
        self.ui.reduceLastButton.hide()
        self.ui.load_live_data_button.hide()

    # Actions defined in Qt Designer
    def file_open_dialog(self):
        """
            Show a dialog to open a new file.
            TODO: consider multiple selection. In this case QuickNXS tries to automatically sort and reduce.
        """

        if self.ui.histogramActive.isChecked():
            filter_=u'Histo Nexus (*histo.nxs);;All (*.*)'
        else:
            filter_=u'Event Nexus (*nxs.h5);;Event Nexus (*event.nxs);;All (*.*)'
        file_path=QtGui.QFileDialog.getOpenFileName(self, u'Open NXS file...',
                                                    directory=self._current_directory,
                                                    filter=filter_)

        if file_path:
            file_dir, file_name = os.path.split(unicode(file_path))
            self.settings.setValue('current_directory', file_dir)

            # Update the list of files
            event_file_list = glob.glob(os.path.join(file_dir, '*event.nxs'))
            h5_file_list = glob.glob(os.path.join(file_dir, '*.nxs.h5'))
            event_file_list.extend(h5_file_list)
            event_file_list.sort()
            event_file_list = [os.path.basename(name) for name in event_file_list]

            current_list=[self.ui.file_list.item(i).text() for i in range(self.ui.file_list.count())]
            if event_file_list != current_list:
                self.ui.file_list.clear()
                for item in event_file_list:
                    listitem = QtGui.QListWidgetItem(item, self.ui.file_list)
                    if item == file_name:
                        self.ui.file_list.setCurrentItem(listitem)
            else:
                try:
                    self.ui.file_list.setCurrentRow(event_file_list.index(file_name))
                except ValueError:
                    logging.error("Could not set file selection %s", file_name)

            self._current_file = file_path
            #self.fileOpen(filenames[0])

    def reload_file(self):
        """
            Reload the file that is currently selected form the list.
        """
        self.fileOpen(self._current_file)


    def plotActiveTab(self): return NotImplemented
    def nextFile(self): return NotImplemented
    def prevFile(self): return NotImplemented
    def setNorm(self): return NotImplemented
    def addRefList(self): return NotImplemented
    def clearRefList(self): return NotImplemented
    def overwriteDirectBeam(self): return NotImplemented
    def normalizeTotalReflection(self): return NotImplemented
    def removeRefList(self): return NotImplemented
    def reduceDatasets(self): return NotImplemented
    def clearOverwrite(self): return NotImplemented
    def loadExtraction(self): return NotImplemented
    def aboutDialog(self): return NotImplemented
    def reductionTableChanged(self): return NotImplemented
    def helpDialog(self): return NotImplemented
    def clearNormList(self): return NotImplemented
    def change_offspec_colorscale(self): return NotImplemented
    def change_gisans_colorscale(self): return NotImplemented
    def open_compare_window(self): return NotImplemented
    def open_advanced_background(self): return NotImplemented
    def clip_offspec_colorscale(self): return NotImplemented
    def fileOpenSumDialog(self): return NotImplemented
    def changeRegionValues(self): return NotImplemented
    def toggleColorbars(self): return NotImplemented
    def overwriteChanged(self): return NotImplemented
    def toggleHide(self): return NotImplemented
    def fileOpenList(self): return NotImplemented
    def visualizePeakfinding(self): return NotImplemented
    def folderModified(self): return NotImplemented
    def replotProjections(self): return NotImplemented
    def openByNumber(self): return NotImplemented
    def cutPoints(self): return NotImplemented
    def autoRef(self): return NotImplemented
    def stripOverlap(self): return NotImplemented
    def changeActiveChannel(self): return NotImplemented
    def open_polarization_window(self): return NotImplemented
    def open_rawdata_dialog(self): return NotImplemented
    def run_ipython(self): return NotImplemented
    def quickReduce(self): return NotImplemented
    def open_nxs_dialog(self): return NotImplemented
    def open_filter_dialog(self): return NotImplemented
    def open_reduction_preview(self): return NotImplemented
    def live_open(self): return NotImplemented
    def open_database_search(self): return NotImplemented
    def exportRawData(self): return NotImplemented


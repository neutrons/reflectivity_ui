#pylint: disable=invalid-name, line-too-long
"""
    Main application window
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
import glob
from PyQt5 import QtGui, QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_main_window

from .data_handling.loader import NexusData
from .configuration import Configuration

class MainWindow(QtWidgets.QMainWindow,
                 reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow):
    """
        Main applicqtion window
    """
    # UI events
    file_loaded_signal = QtCore.pyqtSignal()

    def __init__(self):
        """
            Initialization
        """
        # Base class
        QtWidgets.QMainWindow.__init__(self)

        # Initialize the UI widgets
        self.ui = reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow()
        self.ui.setupUi(self)

        # Application settings
        self.configuration = Configuration()
        self.settings = QtCore.QSettings('.refredm')
        self._current_directory = self.settings.value('current_directory', os.path.expanduser('~'))
        self._current_file = None
        self._active_channel = None

        # Update file list when changes are made        
        self._path_watcher = QtCore.QFileSystemWatcher([self._current_directory], self)
        self._path_watcher.directoryChanged.connect(self.update_file_list)

        # Initialization for specific instrument
        # Retrieve configuration from config and enable/disable features
        self.initialize_instrument()
        self.hide_unsupported()
        
        # UI events
        self.file_loaded_signal.connect(self.update_info)

    def initialize_instrument(self):
        for i in range(1, 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()
        self.ui.selectedChannel0.show()
        self.ui.selectedChannel0.setText(u"none")

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
        self.ui.direct_beam_runs_label.hide()

    # Actions defined in Qt Designer
    def file_open_dialog(self):
        """
            Show a dialog to open a new file.
            TODO: consider multiple selection. In this case QuickNXS tries to automatically sort and reduce.
        """

        if self.ui.histogramActive.isChecked():
            filter_ = u'Histo Nexus (*histo.nxs);;All (*.*)'
        else:
            filter_ = u'Event Nexus (*nxs.h5);;Event Nexus (*event.nxs);;All (*.*)'
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self, u'Open NXS file...',
                                                    directory=self._current_directory,
                                                    filter=filter_)

        if file_path:
            file_dir, file_name = os.path.split(unicode(file_path))
            self.settings.setValue('current_directory', file_dir)

            self._path_watcher.removePath(self._current_directory)
            self._current_directory = file_dir
            self._path_watcher.addPath(self._current_directory)
            self._current_file = file_path
            self._current_file_name = file_name
            self.update_file_list()
            
            # Load the file
            try:
                nexus_data = NexusData(file_path, self.configuration)
                self._data_sets = nexus_data.load()
                self.file_loaded()
            except:
                logging.error("Error loading file: %s", sys.exc_value)
            #self.fileOpen(filenames[0])

    def file_loaded(self):
        current_channel=0
        for i in range(12):
            if getattr(self.ui, 'selectedChannel%i'%i).isChecked():
                current_channel = i

        channels = self._data_sets.keys()
        if current_channel < len(channels):
            self._active_channel = self._data_sets[channels[current_channel]]
        else:
            self._active_channel = self._data_sets[channels[0]]
            self.ui.selectedChannel0.setChecked(True)

        for i, channel in enumerate(channels):
            getattr(self.ui, 'selectedChannel%i'%i).show()
            getattr(self.ui, 'selectedChannel%i'%i).setText(channel)
        for i in range(len(channels), 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()

        #if self.active_channel in self.ref_list_channels:
        #    for i, refli in enumerate(self.reduction_list):
        #        refli=self.recalculateReflectivity(refli)
        #        self.reduction_list[i]=refli
        #self.plotActiveTab()
        #self.initiateProjectionPlot.emit(False)
        #self.initiateReflectivityPlot.emit(False)

        # Update UI
        self.file_loaded_signal.emit()

    def update_info(self):
        '''
        Write file metadata to the labels in the overview tab.
        '''
        d=self._active_channel

        try:
            dangle0=u"%.3f° (%.3f°)"%(float(self.ui.dangle0Overwrite.text()), d.dangle0)
        except ValueError:
            dangle0=u"%.3f°"%(d.dangle0)
        if self.ui.directPixelOverwrite.value()>=0:
            dpix=u"%.1f (%.1f)"%(self.ui.directPixelOverwrite.value(), d.dpix)
        else:
            dpix=u"%.1f"%d.dpix
        self.ui.datasetLambda.setText(u"%.2f (%.2f-%.2f) Å"%(d.lambda_center,
                                                             d.lambda_center-1.5,
                                                             d.lambda_center+1.5))
        self.ui.datasetPCharge.setText(u"%.3e"%d.proton_charge)
        self.ui.datasetTime.setText(u"%i s"%d.total_time)
        self.ui.datasetTotCounts.setText(u"%.4e"%d.total_counts)
        try:
            self.ui.datasetRate.setText(u"%.1f cps"%(d.total_counts/d.total_time))
        except ZeroDivisionError:
            self.ui.datasetRate.setText(u"NaN")
        self.ui.datasetDangle.setText(u"%.3f°"%d.dangle)
        self.ui.datasetDangle0.setText(dangle0)
        self.ui.datasetSangle.setText(u"%.3f°"%d.sangle)
        self.ui.datasetDirectPixel.setText(dpix)
        self.ui.currentChannel.setText('<b>%s</b> (%s)&nbsp;&nbsp;&nbsp;Type: %s&nbsp;&nbsp;&nbsp;Current State: <b>%s</b>'%(
                                       d.number, d.experiment,
                                       d.measurement_type, d.name))

    def update_file_list(self):
        """
            Update the list of data files
        """
        # Update the list of files
        event_file_list = glob.glob(os.path.join(self._current_directory, '*event.nxs'))
        h5_file_list = glob.glob(os.path.join(self._current_directory, '*.nxs.h5'))
        event_file_list.extend(h5_file_list)
        event_file_list.sort()
        event_file_list = [os.path.basename(name) for name in event_file_list]

        current_list=[self.ui.file_list.item(i).text() for i in range(self.ui.file_list.count())]
        if event_file_list != current_list:
            self.ui.file_list.clear()
            for item in event_file_list:
                listitem = QtWidgets.QListWidgetItem(item, self.ui.file_list)
                if item == self._current_file_name:
                    self.ui.file_list.setCurrentItem(listitem)
        else:
            try:
                self.ui.file_list.setCurrentRow(event_file_list.index(self._current_file_name))
            except ValueError:
                logging.error("Could not set file selection: %s", self._current_file_name)
                logging.error(sys.exc_value)

    def reload_file(self):
        """
            Reload the file that is currently selected form the list.
        """
        self.fileOpen(self._current_file)

    def changeActiveChannel(self):
        '''
        The overview and reflectivity channel was changed. This
        recalculates already extracted reflectivities.
        '''
        return self.file_loaded()
    
    def gather_options(self):
        """
            Gather the reduction options.
        """
        pass

    def plotActiveTab(self): return NotImplemented
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


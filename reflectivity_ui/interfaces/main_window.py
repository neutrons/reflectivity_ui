# -*- coding: utf-8 -*-
#pylint: disable=invalid-name, line-too-long, too-many-public-methods, too-many-instance-attributes
"""
    Main application window
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
import glob
import numpy as np
from PyQt5 import QtGui, QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_main_window

from .data_handling.loader import NexusData
from .configuration import Configuration
from . import plotting

class MainWindow(QtWidgets.QMainWindow,
                 reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow):
    """
        Main applicqtion window
    """
    # UI events
    file_loaded_signal = QtCore.pyqtSignal()
    color = None
    _gisansThread = None
    overview_lines = None

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
        self._current_file_name = None
        self._active_channel = None
        self._data_sets = None

        # Update file list when changes are made
        self._path_watcher = QtCore.QFileSystemWatcher([self._current_directory], self)
        self._path_watcher.directoryChanged.connect(self.update_file_list)

        # Initialization for specific instrument
        # Retrieve configuration from config and enable/disable features
        self.initialize_instrument()
        self.hide_unsupported()

        # UI events
        self.file_loaded_signal.connect(self.update_info)
        self.file_loaded_signal.connect(self.update_daslog)
        self.file_loaded_signal.connect(self.plotActiveTab)

    def initialize_instrument(self):
        """
            Initialize instrument according to the instrument
        """
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

    def file_loaded(self):
        """
            Update UI after a file is loaded
        """
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
        """
            Write file metadata to the labels in the overview tab.
        """
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

    def update_daslog(self):
        """
            Write parameters from all file daslogs to the table in the
            daslog tab.
        """
        table=self.ui.daslogTableBox
        table.setRowCount(0)
        table.sortItems(-1)
        table.setColumnCount(len(self._data_sets)+2)
        table.setHorizontalHeaderLabels(['Name']+self._data_sets.keys()+['Unit'])
        for j, key in enumerate(sorted(self._active_channel.logs.keys(), key=lambda s: s.lower())):
            table.insertRow(j)
            table.setItem(j, 0, QtWidgets.QTableWidgetItem(key))
            table.setItem(j, len(self._data_sets)+1,
                          QtWidgets.QTableWidgetItem(self._active_channel.log_units[key]))
            i = 0
            for xs in self._data_sets:
                item=QtWidgets.QTableWidgetItem(u'%g' % self._data_sets[xs].logs[key])
                item.setToolTip(u'MIN: %g   MAX: %g' % (self._data_sets[xs].log_minmax[key]))
                table.setItem(j, i+1, item)
                i += 1
        table.resizeColumnsToContents()

    def plotActiveTab(self):
        '''
            Select the appropriate function to plot all visible images.
        '''
        if self._active_channel is None:
            return
        color=str(self.ui.color_selector.currentText())
        if color!=self.color and self.color is not None:
            self.color=color
            plots=[self.ui.xy_pp, self.ui.xy_mp, self.ui.xy_pm, self.ui.xy_mm,
                   self.ui.xtof_pp, self.ui.xtof_mp, self.ui.xtof_pm, self.ui.xtof_mm,
                   self.ui.xy_overview, self.ui.xtof_overview]
            for plot in plots:
                plot.clear_fig()
        elif self.color is None:
            self.color=color
        if self.ui.plotTab.currentIndex()!=4 and self._gisansThread:
            self._gisansThread.finished.disconnect()
            self._gisansThread.terminate()
            self._gisansThread.wait(100)
            self._gisansThread=None
        if self.ui.plotTab.currentIndex()==0:
            plotting.plot_overview(self)
        if self.ui.plotTab.currentIndex()==1:
            self.plot_xy()
        if self.ui.plotTab.currentIndex()==2:
            self.plot_xtof()
        if self.ui.plotTab.currentIndex()==3:
            self.plot_offspec()
        if self.ui.plotTab.currentIndex()==4:
            self.plot_gisans()

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

    def plot_xy(self):
        '''
        X vs. Y plots for all channels.
        '''
        plots=[self.ui.xy_pp, self.ui.xy_mm, self.ui.xy_pm, self.ui.xy_mp]
        for i in range(len(self.active_data), 4):
            if plots[i].cplot is not None:
                plots[i].clear()
                plots[i].draw()
        imin=1e20
        imax=1e-20
        xynormed=[]
        for dataset in self._data_sets[:4]:
            d=dataset.xydata/dataset.proton_charge
            xynormed.append(d)
            if dataset.total_counts==0:
                continue
            imin=min(imin, d[d>0].min())
            imax=max(imax, d.max())

        if len(xynormed)>1:
            self.ui.frame_xy_mm.show()
            if len(xynormed)==4:
                self.ui.frame_xy_sf.show()
            else:
                self.ui.frame_xy_sf.hide()
        else:
            self.ui.frame_xy_mm.hide()
            self.ui.frame_xy_sf.hide()

        for i, datai in enumerate(xynormed):
            if self.active_data[i].total_counts==0:
                continue
            if self.ui.tthPhi.isChecked():
                plots[i].clear()
                rad_per_pixel=dataset.det_size_x/dataset.dist_sam_det/dataset.xydata.shape[1]
                phi_range=datai.shape[0]*rad_per_pixel*180./np.pi
                tth_range=datai.shape[1]*rad_per_pixel*180./np.pi
                phi0=self.ui.refYPos.value()*rad_per_pixel*180./np.pi
                tth0=(dataset.dangle-dataset.dangle0)-(datai.shape[1]-dataset.dpix)*rad_per_pixel*180./np.pi
                plots[i].imshow(datai, log=self.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, origin='lower',
                                extent=[tth_range+tth0, tth0, phi0, phi0-phi_range])
                plots[i].set_xlabel(u'2$\\Theta{}$ [°]')
                plots[i].set_ylabel(u'$\\phi{}$ [°]')
            else:
                plots[i].imshow(datai, log=self.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, origin='lower')
                plots[i].set_xlabel(u'x [pix]')
                plots[i].set_ylabel(u'y [pix]')
            plots[i].set_title(self.channels[i])
            if plots[i].cplot is not None:
                plots[i].cplot.set_clim([imin, imax])
            if plots[i].cplot is not None and self.ui.show_colorbars.isChecked() and plots[i].cbar is None:
                plots[i].cbar=plots[i].canvas.fig.colorbar(plots[i].cplot)
            plots[i].draw()

    def getNorm(self):
        return None

    def toggleColorbars(self):
        """ Refresh plots because of a color or scale change """
        plots=[self.ui.xy_pp, self.ui.xy_mp, self.ui.xy_pm, self.ui.xy_mm,
               self.ui.xtof_pp, self.ui.xtof_mp, self.ui.xtof_pm, self.ui.xtof_mm,
               self.ui.xy_overview, self.ui.xtof_overview,
               self.ui.offspec_pp, self.ui.offspec_mp, self.ui.offspec_pm, self.ui.offspec_mm]
        for plot in plots:
            plot.clear_fig()
        self.overview_lines=None
        self.plotActiveTab()

    def gather_options(self):
        """
            Gather the reduction options.
        """
        pass

    def setNorm(self): return NotImplemented
    def addRefList(self): return NotImplemented
    def clearRefList(self): return NotImplemented
    def overwriteDirectBeam(self): return NotImplemented
    def normalizeTotalReflection(self): return NotImplemented
    def removeRefList(self): return NotImplemented
    def reduceDatasets(self): return NotImplemented
    def clearOverwrite(self): return NotImplemented
    def loadExtraction(self): return NotImplemented
    def reductionTableChanged(self): return NotImplemented
    def clearNormList(self): return NotImplemented
    def change_offspec_colorscale(self): return NotImplemented
    def change_gisans_colorscale(self): return NotImplemented
    def open_compare_window(self): return NotImplemented
    def open_advanced_background(self): return NotImplemented
    def clip_offspec_colorscale(self): return NotImplemented
    def fileOpenSumDialog(self): return NotImplemented
    def changeRegionValues(self): return NotImplemented
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


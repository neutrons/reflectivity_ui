# -*- coding: utf-8 -*-
"""
    Manage file-related events
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
import glob
import math
from PyQt5 import QtGui, QtCore, QtWidgets


class FileHandler(object):
    def __init__(self, main_window):
        self.ui = main_window.ui
        self.main_window = main_window
        self._data_manager = main_window.data_manager

        # Update file list when changes are made
        self._path_watcher = QtCore.QFileSystemWatcher([self._data_manager.current_directory],
                                                       self.main_window)
        self._path_watcher.directoryChanged.connect(self.update_file_list)

        self.cache_indicator=QtWidgets.QLabel("Cache Size: 0.0MB")
        self.ui.statusbar.addPermanentWidget(self.cache_indicator)
        button=QtWidgets.QPushButton('Empty Cache')
        self.ui.statusbar.addPermanentWidget(button)
        button.pressed.connect(self.empty_cache)
        button.setFlat(True)
        button.setMaximumSize(150, 20)

    def empty_cache(self):
        logging.error("empty cache")
        self._data_manager.clear_cache()
        self.cache_indicator.setText('Cache Size: 0.0MB')

    def open_file(self, file_path, force=False):
        """
            Read a data file
            :param str file_path: file path
            :param bool force: if true, the file will be reloaded
        """
        try:
            self.get_configuration()
            self._data_manager.load(file_path, self.main_window.configuration, force=force)
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

        success = self._data_manager.set_channel(current_channel)
        if not success:
            self.ui.selectedChannel0.setChecked(True)

        channels = self._data_manager.data_sets.keys()
        for i, channel in enumerate(channels):
            getattr(self.ui, 'selectedChannel%i'%i).show()
            getattr(self.ui, 'selectedChannel%i'%i).setText(channel)
        for i in range(len(channels), 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()

        #if self.active_channel in self.ref_list_channels:
        #    for i, refli in enumerate(self.reduction_list):
        #        refli=self.recalculateReflectivity(refli)
        #        self.reduction_list[i]=refli
        #self.initiateReflectivityPlot.emit(False)

        # Update UI
        self.main_window.file_loaded_signal.emit()
        self.main_window.initiate_projection_plot.emit(False)

        self.cache_indicator.setText('Cache Size: %.1fMB'%(self._data_manager.get_cachesize()/1024.**2))

    def update_info(self):
        """
            Write file metadata to the labels in the overview tab.
        """
        d=self._data_manager.active_channel

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

        # Update reduction parameters
        # for lines of the current extraction area
        if d.data_info is not None:
            self.ui.refXPos.setValue(d.data_info.peak_position)
            peak_width = d.data_info.peak_range[1]-d.data_info.peak_range[0]
            self.ui.refXWidth.setValue(peak_width)

            peak_pos = (d.data_info.low_res_range[1]+d.data_info.low_res_range[0])/2.0
            self.ui.refYPos.setValue(peak_pos)
            peak_width = d.data_info.low_res_range[1]-d.data_info.low_res_range[0]
            self.ui.refYWidth.setValue(peak_width)

            bck_pos = (d.data_info.background[1]+d.data_info.background[0])/2.0
            self.ui.bgCenter.setValue(bck_pos)
            bck_width = d.data_info.background[1]-d.data_info.background[0]
            self.ui.bgWidth.setValue(bck_width)

            # TODO: this should update when we change the peak position?
            self.ui.datasetAi.setText(u"%.3f°"%(d.data_info.scattering_angle))
            #self.ui.datasetROI.setText(u"%.4g"%(self.refl.Iraw.sum()))

            self.ui.roi_used_label.setText(u"%s" % d.data_info.use_roi_actual)

        else:
            logging.error("No reduction parameters for this data")

    def update_file_list(self):
        """
            Update the list of data files
        """
        # Update the list of files
        event_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*event.nxs'))
        h5_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*.nxs.h5'))
        event_file_list.extend(h5_file_list)
        event_file_list.sort()
        event_file_list = [os.path.basename(name) for name in event_file_list]

        current_list=[self.ui.file_list.item(i).text() for i in range(self.ui.file_list.count())]
        if event_file_list != current_list:
            self.ui.file_list.clear()
            for item in event_file_list:
                listitem = QtWidgets.QListWidgetItem(item, self.ui.file_list)
                if item == self._data_manager.current_file_name:
                    self.ui.file_list.setCurrentItem(listitem)
        else:
            try:
                self.ui.file_list.setCurrentRow(event_file_list.index(self._data_manager.current_file_name))
            except ValueError:
                logging.error("Could not set file selection: %s", self._data_manager.current_file_name)
                logging.error(sys.exc_value)

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
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self.main_window, u'Open NXS file...',
                                                    directory=self._data_manager.current_directory,
                                                    filter=filter_)

        if file_path:
            file_dir, file_name = os.path.split(unicode(file_path))
            self.main_window.settings.setValue('current_directory', file_dir)

            self._path_watcher.removePath(self._data_manager.current_directory)
            self._data_manager._current_directory = file_dir
            self._data_manager.current_file_name = file_name
            self._path_watcher.addPath(self._data_manager.current_directory)
            self.update_file_list()
            self.open_file(file_path)

    def update_daslog(self):
        """
            Write parameters from all file daslogs to the table in the
            daslog tab.
        """
        table=self.ui.daslogTableBox
        table.setRowCount(0)
        table.sortItems(-1)
        table.setColumnCount(len(self._data_manager.data_sets)+2)
        table.setHorizontalHeaderLabels(['Name']+self._data_manager.data_sets.keys()+['Unit'])
        for j, key in enumerate(sorted(self._data_manager.active_channel.logs.keys(), key=lambda s: s.lower())):
            table.insertRow(j)
            table.setItem(j, 0, QtWidgets.QTableWidgetItem(key))
            table.setItem(j, len(self._data_manager.data_sets)+1,
                          QtWidgets.QTableWidgetItem(self._data_manager.active_channel.log_units[key]))
            i = 0
            for xs in self._data_manager.data_sets:
                item=QtWidgets.QTableWidgetItem(u'%g' % self._data_manager.data_sets[xs].logs[key])
                item.setToolTip(u'MIN: %g   MAX: %g' % (self._data_manager.data_sets[xs].log_minmax[key]))
                table.setItem(j, i+1, item)
                i += 1
        table.resizeColumnsToContents()

    def get_configuration(self):
        """
            Gather the reduction options.
        """
        self.main_window.configuration.tof_bins = self.ui.eventTofBins.value()
        self.main_window.configuration.tof_bin_type = self.ui.eventBinMode.currentIndex()
        self.main_window.configuration.use_roi = self.ui.use_roi_checkbox.isChecked()
        self.main_window.configuration.update_peak_range = self.ui.fit_within_roi_checkbox.isChecked()
        self.main_window.configuration.use_roi_bck = self.ui.use_bck_roi_checkbox.isChecked()

        # Default ranges, using the current values
        x_pos = self.ui.refXPos.value()
        x_width = self.ui.refXWidth.value()
        y_pos = self.ui.refYPos.value()
        y_width = self.ui.refYWidth.value()
        bck_pos = self.ui.bgCenter.value()
        bck_width = self.ui.bgWidth.value()
        
        self.main_window.configuration.forced_peak_roi = [x_pos - x_width/2.0,
                                                          x_pos + x_width/2.0]
        self.main_window.configuration.forced_low_res_roi = [y_pos - y_width/2.0,
                                                             y_pos + y_width/2.0]
        self.main_window.configuration.forced_bck_roi = [bck_pos - bck_width/2.0,
                                                         bck_pos + bck_width/2.0]

        # Use background on each side of the peak
        self.main_window.configuration.use_tight_bck = self.ui.use_side_bck_checkbox.isChecked()
        self.main_window.configuration.bck_offset = self.ui.side_bck_width.value()

        # Other reduction options
        self.main_window.configuration.subtract_background = self.ui.bgActive.isChecked()
        self.main_window.configuration.scaling_factor = self.ui.refScale.value()
        self.main_window.configuration.cut_first_n_points = self.ui.rangeStart.value()
        self.main_window.configuration.cut_last_n_points = self.ui.rangeEnd.value()

        # Make the changes persistent
        self.main_window.configuration.to_q_settings(self.main_window.settings)

    def populate_from_configuration(self):
        """
            Set reduction options in UI, usually after loading
            a reduced data set.
        """
        self.ui.eventTofBins.setValue(self.main_window.configuration.tof_bins)
        self.ui.eventBinMode.setCurrentIndex(self.main_window.configuration.tof_bin_type)
        self.ui.use_roi_checkbox.setChecked(self.main_window.configuration.use_roi)
        self.ui.fit_within_roi_checkbox.setChecked(self.main_window.configuration.update_peak_range)
        self.ui.use_bck_roi_checkbox.setChecked(self.main_window.configuration.use_roi_bck)

        # Default ranges, using the current values
        x_pos = (self.main_window.configuration.forced_peak_roi[1] \
                 + self.main_window.configuration.forced_peak_roi[0]) / 2.0
        x_width = (self.main_window.configuration.forced_peak_roi[1] \
                   - self.main_window.configuration.forced_peak_roi[0])

        y_pos = (self.main_window.configuration.forced_low_res_roi[1] \
                 + self.main_window.configuration.forced_low_res_roi[0]) / 2.0
        y_width = (self.main_window.configuration.forced_low_res_roi[1] \
                   - self.main_window.configuration.forced_low_res_roi[0])

        bck_pos = (self.main_window.configuration.forced_bck_roi[1] \
                   + self.main_window.configuration.forced_bck_roi[0]) / 2.0
        bck_width = (self.main_window.configuration.forced_bck_roi[1] \
                     - self.main_window.configuration.forced_bck_roi[0])

        self.ui.refXPos.setValue(x_pos)
        self.ui.refXWidth.setValue(x_width)
        self.ui.refYPos.setValue(y_pos)
        self.ui.refYWidth.setValue(y_width)
        self.ui.bgCenter.setValue(bck_pos)
        self.ui.bgWidth.setValue(bck_width)

        # Use background on each side of the peak
        self.ui.use_side_bck_checkbox.setChecked(self.main_window.configuration.use_tight_bck)
        self.ui.side_bck_width.setValue(self.main_window.configuration.bck_offset)

        # Subtract background
        self.ui.bgActive.setChecked(self.main_window.configuration.subtract_background)
        # Scaling factor
        self.ui.refScale.setValue(math.log10(self.main_window.configuration.scaling_factor))
        # Cut first and last points
        self.ui.rangeStart.setValue(self.main_window.configuration.cut_first_n_points)
        self.ui.rangeEnd.setValue(self.main_window.configuration.cut_last_n_points)

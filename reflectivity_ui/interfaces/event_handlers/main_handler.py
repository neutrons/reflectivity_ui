# -*- coding: utf-8 -*-
#pylint: disable=invalid-name, line-too-long, too-many-public-methods, too-many-instance-attributes, wrong-import-order, bare-except, protected-access, too-many-arguments, too-many-statements
"""
    Manage file-related and UI events
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
import glob
import math
import time
from PyQt5 import QtGui, QtCore, QtWidgets

from ..configuration import Configuration
from .progress_reporter import ProgressReporter


class MainHandler(object):
    """
        Event handler for the main application window.
    """
    def __init__(self, main_window):
        self.ui = main_window.ui
        self.main_window = main_window
        self._data_manager = main_window.data_manager

        # Update file list when changes are made
        self._path_watcher = QtCore.QFileSystemWatcher([self._data_manager.current_directory],
                                                       self.main_window)
        self._path_watcher.directoryChanged.connect(self.update_file_list)

        self.cache_indicator = QtWidgets.QLabel("Files loaded: 0")
        self.cache_indicator.setMargin(5)
        self.cache_indicator.setSizePolicy(QtWidgets.QSizePolicy.Fixed,
                                           QtWidgets.QSizePolicy.Preferred)
        self.cache_indicator.setMinimumWidth(110)
        self.ui.statusbar.addPermanentWidget(self.cache_indicator)
        button = QtWidgets.QPushButton('Empty Cache')
        self.ui.statusbar.addPermanentWidget(button)
        button.pressed.connect(self.empty_cache)
        button.setFlat(True)
        button.setMaximumSize(150, 20)

        # Create progress bar in statusbar
        self.progress_bar = QtWidgets.QProgressBar(self.ui.statusbar)
        self.progress_bar.setMinimumSize(20, 14)
        self.progress_bar.setMaximumSize(140, 100)
        self.ui.statusbar.addPermanentWidget(self.progress_bar)

        self.status_message = QtWidgets.QLabel("")
        self.status_message.setMinimumWidth(1000)
        self.status_message.setMargin(5)
        self.ui.statusbar.insertWidget(0, self.status_message)

    def new_progress_reporter(self):
        """ Return a progress reporter """
        return ProgressReporter(progress_bar=self.progress_bar, status_bar=self.status_message)

    def empty_cache(self):
        """
            Empty the data cache
        """
        self._data_manager.clear_cache()
        self.cache_indicator.setText("Files loaded: 0")

    def open_file(self, file_path, force=False, silent=False):
        """
            Read a data file
            :param str file_path: file path
            :param bool force: if true, the file will be reloaded
            :param bool silent: if true, the UI will not be updated
        """
        print('[DEBUG] open_file is called for {}'.format(file_path))
        if not os.path.isfile(file_path):
            self.report_message("File does not exist",
                                detailed_message="The following file does not exist:\n  %s" % file_path,
                                pop_up=True, is_error=True)
            return
        t_0 = time.time()
        self.main_window.auto_change_active = True
        try:
            self.report_message("Loading file %s" % file_path)
            prog = ProgressReporter(progress_bar=self.progress_bar, status_bar=self.status_message)
            configuration = self.get_configuration()
            self._data_manager.load(file_path, configuration, force=force, progress=prog)
            self.report_message("Loaded file %s" % self._data_manager.current_file_name)
        except RuntimeError as run_err:
            # FIXME - need to find out what kind of error it could have
            self.report_message("Error loading file {} due to {}".format(self._data_manager.current_file_name, run_err),
                                detailed_message=str(sys.exc_value), pop_up=False, is_error=True)
        # FIXME 63 - Disable general exception for future
	# except Exception as e:
        #     print('Unhandled general exception {}'.format(e))
        #     raise e
       
        if not silent:
            self.file_loaded()
        self.main_window.auto_change_active = False
        logging.info("DONE: %s sec", time.time()-t_0)

    def load_merge_files(self, file_paths, force=False, silent=False):
        """Load and merge multiple Nexus file

        Event handling workflow
        - reflectivity_ui/interfaces/main_window.py", line 136, in files_open_dialog
          self.file_handler.files_open_dialog()
        - reflectivity_ui/interfaces/event_handlers/main_handler.py", line 515, in files_open_dialog
          self.load_merge_files(file_paths)

        Examples:
        - file_path: [u'/SNS/REF_M/IPTS-25531/nexus/REF_M_38189.nxs.h5',
                      u'/SNS/REF_M/IPTS-25531/nexus/REF_M_38189.nxs.h5']
        -

        Parameters
        ----------
        file_paths: ~list
            list of Nexus file path with full path
        force: bool
            flag to force
        silent: bool
            flag to be silent.  If not explicitly call method file_loaded()

        Returns
        -------
        None

        """
        # TODO 63: TASK 1 (in progress): clean, document and set fake
        # Verify that all files selected shall exist
        for file_path in file_paths:
            if not os.path.isfile(file_path):
                self.report_message("File does not exist",
                                    detailed_message="The following file does not exist:\n  %s" % file_path,
                                    pop_up=True, is_error=True)
                return

        # Start timing
        t_0 = time.time()
        self.main_window.auto_change_active = True
        try:
            self.report_message("Loading files {}".format(file_paths))
            prog = ProgressReporter(progress_bar=self.progress_bar, status_bar=self.status_message)
            configuration = self.get_configuration()
            # FIXME 63 (clean and doc) data manage of of type
            #  <class 'reflectivity_ui.interfaces.data_manager.DataManager'>
            self._data_manager.load_merge(file_paths, configuration, force=force, progress=prog)
            self.report_message("Loaded file %s" % self._data_manager.current_file_name)
        except RuntimeError as e:
            # Loading error: report and return
            self.report_message("Error loading file %s" % self._data_manager.current_file_name,
                                detailed_message=str(sys.exc_value), pop_up=False, is_error=True)
            return
        except Exception as e:
            print('General exception {} is not handled well'.format(e))
            raise e

        # TODO 66: Task 66/67 takes over from here --------------------------
        # response for file loaded
        if not silent:
            self.merged_files_loaded()

        # set flag
        self.main_window.auto_change_active = False
        logging.info("DONE: %s sec", time.time()-t_0)

    def merged_files_loaded(self):
        """
            Update UI after a file is loaded

            communication is via self._data_manager
        """
        # TODO 66 - Implement ASAP by referring to method 'file_loaded()'
        #  66 How to make it work with merged data???
        self.main_window.auto_change_active = True
        current_channel = 0

        # TODO FIXME Is 12 a magic number?
        # Get selected channel from UI
        # FIXME - observe the UI!
        for i in range(12):
            if getattr(self.ui, 'selectedChannel%i' % i).isChecked():
                current_channel = i
        print('[DEBUG 66] select channel {}'.format(current_channel))

        # Set channel to data manager
        success = self._data_manager.set_channel(current_channel)
        if not success:
            self.ui.selectedChannel0.setChecked(True)

        print('[DEBUG] data manager type: {}'.format(type(self._data_manager)))
        print('[DEBUG] data manager data set keys: {}'.format(list(self._data_manager.data_sets.keys())))
        channels = self._data_manager.data_sets.keys()
        for i, channel in enumerate(channels):
            getattr(self.ui, 'selectedChannel%i'%i).show()
            good_label = channel.replace('_', '-')
            if not good_label == self._data_manager.data_sets[channel].cross_section_label:
                good_label = "%s: %s" % (good_label, self._data_manager.data_sets[channel].cross_section_label)
            getattr(self.ui, 'selectedChannel%i'%i).setText(good_label)
        for i in range(len(channels), 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()
        self.main_window.auto_change_active = False

        # TODO #66
        #    --->  reflectivity_ui/interfaces/event_handlers/main_handler.py:: update_info(self):
        self.main_window.file_loaded_signal.emit()
        self.main_window.initiate_reflectivity_plot.emit(False)
        self.main_window.initiate_projection_plot.emit(False)

        self.cache_indicator.setText('Files loaded: %s' % (self._data_manager.get_cachesize()))

    def file_loaded(self):
        """
            Update UI after a file is loaded
        """
        self.main_window.auto_change_active = True
        current_channel = 0
        for i in range(12):
            if getattr(self.ui, 'selectedChannel%i'%i).isChecked():
                current_channel = i

        success = self._data_manager.set_channel(current_channel)
        if not success:
            self.ui.selectedChannel0.setChecked(True)

        channels = self._data_manager.data_sets.keys()
        for i, channel in enumerate(channels):
            getattr(self.ui, 'selectedChannel%i'%i).show()
            good_label = channel.replace('_', '-')
            if not good_label == self._data_manager.data_sets[channel].cross_section_label:
                good_label = "%s: %s" % (good_label, self._data_manager.data_sets[channel].cross_section_label)
            getattr(self.ui, 'selectedChannel%i'%i).setText(good_label)
        for i in range(len(channels), 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()
        self.main_window.auto_change_active = False

        self.main_window.file_loaded_signal.emit()
        self.main_window.initiate_reflectivity_plot.emit(False)
        self.main_window.initiate_projection_plot.emit(False)

        self.cache_indicator.setText('Files loaded: %s' % (self._data_manager.get_cachesize()))

    def check_files_to_merge(self, file_paths):
        """Check whether these files can be merged

        Parameters
        ----------
        file_paths: ~list
            List of Nexus files (full path)

        Returns
        -------
        str, None
           Error message

        """
        message = self._data_manager.check_files_for_merging(file_paths)

        return message

    def update_tables(self):
        """
            Update a data set that may be in the reduction table or the
            direct beam table.
        """
        # Update the reduction table if this data set is in it
        idx = self._data_manager.find_active_data_id()
        if idx is not None:
            self.update_reduction_table(idx, self._data_manager.active_channel)

        # Update the direct beam table if this data set is in it
        idx = self._data_manager.find_active_direct_beam_id()
        if idx is not None:
            self.update_direct_beam_table(idx, self._data_manager.active_channel)

    def update_calculated_data(self):
        """
            Update the calculated entries in the overview tab.
            We should call this after the peak ranges change, or
            after a change is made that will affect the displayed results.
        """
        d = self._data_manager.active_channel
        self.ui.datasetAi.setText(u"%.3f°"%(d.scattering_angle))

        wl_min, wl_max = d.wavelength_range
        self.ui.datasetLambda.setText(u"%.2f (%.2f-%.2f) Å"%(d.lambda_center,
                                                             wl_min, wl_max))

        # DIRPIX and DANGLE0 overwrite
        if self.ui.set_dangle0_checkbox.isChecked():
            dangle0 = u"%.3f° (%.3f°)" % (float(self.ui.dangle0Overwrite.text()), d._angle_offset)
        else:
            dangle0 = u"%.3f°"%(d.angle_offset)
        self.ui.datasetDangle0.setText(dangle0)

        if self.ui.set_dirpix_checkbox.isChecked():
            dpix = u"%.1f (%.1f)" % (float(self.ui.directPixelOverwrite.value()), d._direct_pixel)
        else:
            dpix = u"%.1f"%d.direct_pixel
        self.ui.datasetDirectPixel.setText(dpix)

        if d.configuration.normalization is not None:
            self.ui.matched_direct_beam_label.setText(u"%s" % d.configuration.normalization)
        else:
            self.ui.matched_direct_beam_label.setText(u"None")

    def update_info(self):
        """
            Update metadata shown in the overview tab.
        """
        self.main_window.auto_change_active = True
        d = self._data_manager.active_channel

        if d is None:
            # TODO 66 - find out what if d is None
            return
            # raise RuntimeError('active channel d is None')
        else:
            print('active channel d is of type {}'.format(type(d)))

        self.populate_from_configuration(d.configuration)
        self.main_window.initiate_projection_plot.emit(False)
        QtWidgets.QApplication.instance().processEvents()

        if self.ui.set_dangle0_checkbox.isChecked():
            dangle0 = u"%.3f° (%.3f°)" % (float(self.ui.dangle0Overwrite.text()), d.angle_offset)
        else:
            dangle0 = u"%.3f°"%(d.angle_offset)

        if self.ui.set_dirpix_checkbox.isChecked():
            dpix = u"%.1f (%.1f)" % (float(self.ui.directPixelOverwrite.value()), d.direct_pixel)
        else:
            dpix = u"%.1f"%d.direct_pixel

        wl_min, wl_max = d.wavelength_range
        self.ui.datasetLambda.setText(u"%.2f (%.2f-%.2f) Å"%(d.lambda_center,
                                                             wl_min, wl_max))
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
        self.ui.currentChannel.setText('<b>%s</b> (%s)&nbsp;&nbsp;&nbsp;Type: %s&nbsp;&nbsp;&nbsp;Current State: '
                                       '<b>%s</b>'%(d.number, d.experiment, d.measurement_type, d.name))

        # Update direct beam indicator
        if d.is_direct_beam:
            self.ui.is_direct_beam_label.setText(u"Direct beam")
        else:
            self.ui.is_direct_beam_label.setText(u"")

        # Update the calculated data
        self.update_calculated_data()

        self.ui.roi_used_label.setText(u"%s" % d.use_roi_actual)
        self.ui.roi_peak_label.setText(u"%s" % str(d.meta_data_roi_peak))
        self.ui.roi_bck_label.setText(u"%s" % str(d.meta_data_roi_bck))

        # Update reduction tables
        self.update_tables()

        self.active_data_changed()

        self.main_window.auto_change_active = False

    def update_file_list(self, file_path=None, merge_mode=False):
        """Update the list of data files

        Parameters
        ----------
        file_path: str, ~list
            Full nexus file path or a list of nexut file path
        merge_mode: bool
            Flag if it is in merge mode

        Returns
        -------
        None

        """
        # TODO 63 [in progress] #63 - Implement if merge_mode is True

        # set the UI mode
        self.main_window.auto_change_active = True

        # TODO 63 (new): merge mode
        if merge_mode:
            # check file
            if isinstance(file_path, list) is False:
                raise RuntimeError('Must be list')

            # set path
            file_paths = file_path
            file_path = file_paths[0]
        else:
            file_paths = None

        if file_path is not None and not file_path == self._data_manager.current_directory:
            print('[DEBUG] File path is not None and data manager current directory = {}'.format(file_path))
            if os.path.isdir(file_path):
                file_dir = file_path
            else:
                file_dir, file_name = os.path.split(unicode(file_path))
                self._data_manager.current_file_name = file_name
            self.main_window.settings.setValue('current_directory', file_dir)
            self._path_watcher.removePath(self._data_manager.current_directory)
            self._data_manager.current_directory = file_dir
            self._path_watcher.addPath(self._data_manager.current_directory)
        else:
            # do nothing
            print('[DEBUG] {} is either None or not equal to data manager current directory "{}"'
                  ''.format(file_path, self._data_manager.current_directory))

        # Update the list of files
        event_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*event.nxs'))
        h5_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*.nxs.h5'))
        event_file_list.extend(h5_file_list)
        event_file_list.sort()
        event_file_list = [os.path.basename(name) for name in event_file_list]

        current_list = [self.ui.file_list.item(i).text() for i in range(self.ui.file_list.count())]
        if event_file_list != current_list:
            # Reset ui.file_list
            print('[DEBUG UI] current list = {} Not equal to event lsit'.format(current_list))
            self.ui.file_list.clear()

            # TODO 63: combined item (in progress)
            if merge_mode:
                # combine file names for a new entry: as f1+f2+
                merged_item = ''
                for fi, item in enumerate(file_paths):
                    if fi > 0:
                        merged_item += '+'
                    merged_item += item.split('.')[0]
                # add to UI file list
                QtWidgets.QListWidgetItem(merged_item, self.ui.file_list)
                # FIXME 63: does this work?
                self.ui.file_list.setCurrentItem(merged_item)

            # add rest of the files
            print('[DEBUG 63-64] data manager current file name: {}'.format(self._data_manager.current_file_name))
            for item in event_file_list:
                # create and add QListWidgetItem from self.ui.file_list by 
                listitem = QtWidgets.QListWidgetItem(item, self.ui.file_list)
                # current_file_name: focus the file list table
                if item == self._data_manager.current_file_name:
                    self.ui.file_list.setCurrentItem(listitem)
        else:
            # Focus ui.file_list to current file name
            try:
                self.ui.file_list.setCurrentRow(event_file_list.index(self._data_manager.current_file_name))
            except ValueError:
                self.report_message("Could not set file selection: %s" % self._data_manager.current_file_name,
                                    detailed_message=str(sys.exc_value), pop_up=False, is_error=True)

        # reset the flag
        self.main_window.auto_change_active = False

    def automated_file_selection(self):
        """
            Go through the files in the current in order of run numbers, and
            load files until the incident angle is no longer increasing.
        """
        self.main_window.auto_change_active = True
        # Update the list of files
        event_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*event.nxs'))
        h5_file_list = glob.glob(os.path.join(self._data_manager.current_directory, '*.nxs.h5'))
        event_file_list.extend(h5_file_list)
        event_file_list.sort()
        event_file_list = [os.path.basename(name) for name in event_file_list]

        current_file_found = False
        n_count = 0
        logging.error("Current file: %s", self._data_manager.current_file_name)

        q_current = self._data_manager.extract_meta_data().mid_q

        # Add the current data set to the reduction table
        # Do nothing if the data is incompatible
        is_direct_beam = self._data_manager.active_channel.is_direct_beam
        if is_direct_beam:
            if not self.add_direct_beam():
                return
        else:
            if not self.add_reflectivity():
                return

        for f in event_file_list:
            file_path = str(os.path.join(self._data_manager.current_directory, f))
            if current_file_found and n_count < 10:
                n_count += 1
                meta_data = self._data_manager.extract_meta_data(file_path)

                if q_current <= meta_data.mid_q and is_direct_beam == meta_data.is_direct_beam:
                    q_current = meta_data.mid_q
                    self.open_file(file_path, silent=True)
                    d = self._data_manager.active_channel
                    # If we find data of another type, stop here
                    if not is_direct_beam == self._data_manager.active_channel.is_direct_beam:
                        break
                    self.main_window.auto_change_active = True
                    self.populate_from_configuration(d.configuration)
                    if self._data_manager.active_channel.is_direct_beam:
                        self.add_direct_beam()
                    else:
                        self.add_reflectivity()

            if f == self._data_manager.current_file_name:
                current_file_found = True

        # At the very end, update the UI and plot reflectivity
        if n_count > 0:
            self.main_window.auto_change_active = True
            self.file_loaded()
        self.main_window.auto_change_active = False

    def open_reduced_file_dialog(self):
        """
            Open a reduced file and all the data files needed to reproduce it.
        """
        # Open file dialog
        filter_ = u'QuickNXS files (*.dat);;All (*.*)'
        output_dir = self.main_window.settings.value('output_directory', os.path.expanduser('~'))
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self.main_window, u'Open reduced file...',
                                                             directory=output_dir,
                                                             filter=filter_)

        t_0 = time.time()
        if file_path:
            # Clear the reduction list first so that we don't create problems later
            self.clear_direct_beams()
            self.clear_reflectivity()
            configuration = self.get_configuration()
            prog = self.new_progress_reporter()
            self._data_manager.load_data_from_reduced_file(file_path, configuration=configuration,
                                                           progress=prog)

            # Update output directory
            file_dir, _ = os.path.split(unicode(file_path))
            self.main_window.settings.setValue('output_directory', file_dir)

            self.main_window.auto_change_active = True

            self.ui.normalizeTable.setRowCount(len(self._data_manager.direct_beam_list))
            for idx, _ in enumerate(self._data_manager.direct_beam_list):
                self._data_manager.set_active_data_from_direct_beam_list(idx)
                self.update_direct_beam_table(idx, self._data_manager.active_channel)
            self.ui.reductionTable.setRowCount(len(self._data_manager.reduction_list))
            for idx, _ in enumerate(self._data_manager.reduction_list):
                self._data_manager.set_active_data_from_reduction_list(idx)
                self.update_reduction_table(idx, self._data_manager.active_channel)

            direct_beam_ids = [str(r.number) for r in self._data_manager.direct_beam_list]
            self.ui.normalization_list_label.setText(u", ".join(direct_beam_ids))

            self.file_loaded()

            if self._data_manager.active_channel is not None:
                self.populate_from_configuration(self._data_manager.active_channel.configuration)
                self.update_file_list(self._data_manager.current_file)
            self.main_window.auto_change_active = False

            logging.info("UI updated: %s", time.time()-t_0)

    # Actions defined in Qt Designer
    def file_open_dialog(self):
        """
            Show a dialog to open a new file.
            TODO: consider multiple selection. In this case QuickNXS tries to automatically sort and reduce.
        """
        if self.ui.histogramActive.isChecked():
            filter_ = u'All (*.*);;histo.nxs (*histo.nxs)'
        else:
            filter_ = u'All (*.*);;nxs.h5 (*nxs.h5);;event.nxs (*event.nxs)'
        file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self.main_window, u'Open NXS file...',
                                                             directory=self._data_manager.current_directory,
                                                             filter=filter_)

        if file_path:
            self.update_file_list(file_path)
            self.open_file(file_path)

    # Actions defined in Qt Designer
    def files_open_dialog(self):
        """
            Show a dialog to open a new file.
            TODO: consider multiple selection. In this case QuickNXS tries to automatically sort and reduce.
        """
        def select_files():
            """launch a dialog to select multiple files

            Returns
            -------
            ~list
                list of full path to Nexus files selected

            """
            # TODO 63 - Implement the dialog to select multiple files
            # file_path, _ = QtWidgets.QFileDialog.getOpenFileName(self.main_window, u'Open NXS file...',
            #                                                      directory=self._data_manager.current_directory,
            #                                                      filter=filter_)

            # FIXME 63 - using fake data now
            selected_file_paths = ['/SNS/REF_M/IPTS-25531/nexus/REF_M_38189.nxs.h5',
                                   '/SNS/REF_M/IPTS-25531/nexus/REF_M_38189.nxs.h5']
            return selected_file_paths


        # Set file filters
        if self.ui.histogramActive.isChecked():
            filter_ = u'All (*.*);;histo.nxs (*histo.nxs)'
        else:
            filter_ = u'All (*.*);;nxs.h5 (*nxs.h5);;event.nxs (*event.nxs)'
        # FIXME TODO - Replace by multiple files selector
        file_paths = select_files()

        # user cancel operation
        if len(file_paths) == 0:
            return

        # Process files that are selected
        # check whether files can be merged without further notice
        message = self.check_files_to_merge(file_paths)
        # need user's permission
        if message or len(message) > 0:
            # TODO # 65 - Implement method 'ask_user_permission'
            user_say_go = self.ask_user_permission(message)
            if not user_say_go:
                return
        # udpate file list
        self.update_file_list(file_paths, merge_mode=True)
        # merge and load file
        self.load_merge_files(file_paths)

    def ask_user_permission(self, message):
        """Ask user's permission to proceed or return
        """
        # TODO 65 - Implement
        print('[DEBUG] Show message: "{}" and ask user to proceed or not.'.format(message))

        return True

    def open_run_number(self, number=None):
        """
            Open a data file by typing a run number
        """
        self.main_window.auto_change_active = True
        if number is None:
            number = self.ui.numberSearchEntry.text()
            # TODO FIXME number can be numbers
        QtWidgets.QApplication.instance().processEvents()

        # Look for new-style nexus file name
        configuration = self.get_configuration()
        search_string = configuration.instrument.file_search_template % number

        file_list = glob.glob(search_string+'.nxs.h5')
        # Look for old-style nexus file name
        if not file_list:
            search_string = configuration.instrument.legacy_search_template % number
            file_list = glob.glob(search_string+'_event.nxs')
        self.ui.numberSearchEntry.setText('')

        success = False
        if file_list > 0:
            self.update_file_list(file_list[0])
            self.open_file(os.path.abspath(file_list[0]))
            success = True
        else:
            self.report_message("Could not locate file %s" % number, pop_up=True)

        self.main_window.auto_change_active = False
        return success

    def update_daslog(self):
        """
            Write parameters from all file daslogs to the table in the
            daslog tab.
        """
        table = self.ui.daslogTableBox
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
                item = QtWidgets.QTableWidgetItem(u'%g' % self._data_manager.data_sets[xs].logs[key])
                item.setToolTip(u'MIN: %g   MAX: %g' % (self._data_manager.data_sets[xs].log_minmax[key]))
                table.setItem(j, i+1, item)
                i += 1
        table.resizeColumnsToContents()

    def add_reflectivity(self, silent=False):
        """
            Collect information about the current extraction settings and store them
            in the list of reduction items.

            Returns true if everything is ok, false otherwise.
        """
        # Update the configuration according to current parameters
        # Note that when a data set is first loaded, the peaks may have a different
        # range for each cross-section. If the option to use a common set of ranges
        # was turned on, we pick the ranges from the currently active cross-section
        # and apply then to all cross-sections.
        if self.ui.action_use_common_ranges.isChecked():
            config = self.get_configuration()
            self._data_manager.update_configuration(configuration=config, active_only=False)

        # Verify that the new data is consistent with existing data in the table
        if not self._data_manager.add_active_to_reduction():
            if not silent:
                self.report_message("Data incompatible or already in the list.", pop_up=True)
            return False
        self.main_window.auto_change_active = True

        # Update the reduction and direct beam tables
        idx = self._data_manager.find_data_in_reduction_list(self._data_manager._nexus_data)
        self.ui.reductionTable.insertRow(idx)
        self.update_tables()

        self.main_window.initiate_reflectivity_plot.emit(True)
        self.main_window.update_specular_viewer.emit()
        self.main_window.auto_change_active = False
        return True

    def update_reduction_table(self, idx, d):
        """
            Update the reduction tale
        """
        self.main_window.auto_change_active = True
        item = QtWidgets.QTableWidgetItem(str(d.number))
        if d == self._data_manager.active_channel:
            item.setBackground(QtGui.QColor(246, 213, 16))
        else:
            item.setBackground(QtGui.QColor(255, 255, 255))
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
        self.ui.reductionTable.setItem(idx, 0, item)
        self.ui.reductionTable.setItem(idx, 1,
                                       QtWidgets.QTableWidgetItem("%.4f"%(d.configuration.scaling_factor)))
        self.ui.reductionTable.setItem(idx, 2,
                                       QtWidgets.QTableWidgetItem(str(d.configuration.cut_first_n_points)))
        self.ui.reductionTable.setItem(idx, 3,
                                       QtWidgets.QTableWidgetItem(str(d.configuration.cut_last_n_points)))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.peak_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.reductionTable.setItem(idx, 4, item)
        self.ui.reductionTable.setItem(idx, 5,
                                       QtWidgets.QTableWidgetItem(str(d.configuration.peak_width)))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.low_res_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.reductionTable.setItem(idx, 6, item)
        self.ui.reductionTable.setItem(idx, 7,
                                       QtWidgets.QTableWidgetItem(str(d.configuration.low_res_width)))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.bck_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.reductionTable.setItem(idx, 8, item)
        self.ui.reductionTable.setItem(idx, 9,
                                       QtWidgets.QTableWidgetItem(str(d.configuration.bck_width)))
        self.ui.reductionTable.setItem(idx, 10,
                                       QtWidgets.QTableWidgetItem(str(d.direct_pixel)))
        self.ui.reductionTable.setItem(idx, 11,
                                       QtWidgets.QTableWidgetItem("%.4f"%d.scattering_angle))
        norma = 'none'
        if d.configuration.normalization is not None:
            norma = d.configuration.normalization
        self.ui.reductionTable.setItem(idx, 12,
                                       QtWidgets.QTableWidgetItem(str(norma)))
        self.main_window.auto_change_active = False

    def clear_reflectivity(self):
        """
            Remove all items from the reduction list.
        """
        self._data_manager.reduction_list = []
        self.ui.reductionTable.setRowCount(0)
        self.main_window.initiate_reflectivity_plot.emit(False)

    def clear_direct_beams(self):
        """
            Remove all items from the direct beam list.
        """
        self._data_manager.clear_direct_beam_list()
        self.ui.normalizeTable.setRowCount(0)
        self.ui.normalization_list_label.setText(u"None")
        self.main_window.initiate_reflectivity_plot.emit(False)

    def remove_reflectivity(self):
        """
            Remove one item from the reduction list.
        """
        index = self.ui.reductionTable.currentRow()
        if index < 0:
            return
        self._data_manager.reduction_list.pop(index)
        self.ui.reductionTable.removeRow(index)
        self.main_window.initiate_reflectivity_plot.emit(False)

    def remove_direct_beam(self):
        """
            Remove one item from the direct beam list.
        """
        index = self.ui.normalizeTable.currentRow()
        if index < 0:
            return
        self._data_manager.direct_beam_list.pop(index)
        self.ui.normalizeTable.removeRow(index)
        self.main_window.initiate_reflectivity_plot.emit(False)

    def reduction_table_changed(self, item):
        '''
            Perform action upon change in data reduction list.
        '''
        if self.main_window.auto_change_active:
            return

        entry = item.row()
        column = item.column()

        refl = self._data_manager.reduction_list[entry]

        #TODO: If we changed the normalization run, make sure it's in the list
        # of direct beams we know about.

        keys = ['number', 'scaling_factor', 'cut_first_n_points', 'cut_last_n_points',
                'peak_position', 'peak_width', 'low_res_position', 'low_res_width',
                'bck_position', 'bck_width', 'direct_pixel', 'scattering_angle', 'normalization']

        # Update settings from selected option
        if column in [1, 4, 5, 6, 7, 8, 9, 10]:
            refl.set_parameter(keys[column], float(item.text()))
        elif column in [2, 3]:
            refl.set_parameter(keys[column], int(item.text()))
        elif column == 12:
            try:
                refl.set_parameter(keys[column], item.text())
            except:
                refl.set_parameter(keys[column], None)
                item.setText("none")

        # Update calculated data
        refl.update_calculated_values()

        # If the changed data set is the active data, also change the UI
        if self._data_manager.is_active(refl):
            self.main_window.auto_change_active = True
            self.update_info()
            self.main_window.auto_change_active = False

        # Update the direct beam table if this data set is in it
        idx = self._data_manager.find_data_in_direct_beam_list(refl)
        if idx is not None:
            channels = refl.cross_sections.keys()
            self.update_direct_beam_table(idx, refl.cross_sections[channels[0]])

        # Only recalculate if we need to, otherwise just replot
        if not column in [1, 2, 3]:
            try:
                self._data_manager.calculate_reflectivity(nexus_data=refl)
            except:
                self.report_message("Could not compute reflectivity for %s" % self._data_manager.current_file_name,
                                    detailed_message=str(sys.exc_value), pop_up=False, is_error=False)

        self.main_window.initiate_reflectivity_plot.emit(True)
        self.main_window.update_specular_viewer.emit()

    def add_direct_beam(self, silent=False):
        """
            Add / remove dataset to the available normalizations or clear the normalization list.
        """
        # Update all cross-section parameters as needed.
        if self.ui.action_use_common_ranges.isChecked():
            config = self.get_configuration()
            self._data_manager.update_configuration(configuration=config, active_only=False)

        # Verify that the new data is consistent with existing data in the table
        if not self._data_manager.add_active_to_normalization():
            if not silent:
                self.report_message("Data incompatible or already in the list.", pop_up=True)
            return False

        self.ui.normalizeTable.setRowCount(len(self._data_manager.direct_beam_list))
        self.update_tables()

        direct_beam_ids = [str(r.number) for r in self._data_manager.direct_beam_list]
        self.ui.normalization_list_label.setText(u", ".join(direct_beam_ids))

        self.main_window.initiate_reflectivity_plot.emit(False)
        return True

    def update_direct_beam_table(self, idx, d):
        """
            Update a direct beam table entry
            :param int idx: row index
            :param CrossSectionData d: data object
        """
        self.main_window.auto_change_active = True
        item = QtWidgets.QTableWidgetItem(str(d.number))
        item.setFlags(item.flags() & ~QtCore.Qt.ItemIsEditable)
        if d == self._data_manager.active_channel:
            item.setBackground(QtGui.QColor(246, 213, 16))
        else:
            item.setBackground(QtGui.QColor(255, 255, 255))

        self.ui.normalizeTable.setItem(idx, 0, QtWidgets.QTableWidgetItem(item))
        wl = u"%s - %s" % (d.wavelength[0], d.wavelength[-1])
        self.ui.normalizeTable.setItem(idx, 7, QtWidgets.QTableWidgetItem(wl))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.peak_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.normalizeTable.setItem(idx, 1, QtWidgets.QTableWidgetItem(item))
        self.ui.normalizeTable.setItem(idx, 2, QtWidgets.QTableWidgetItem(str(d.configuration.peak_width)))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.low_res_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.normalizeTable.setItem(idx, 3, QtWidgets.QTableWidgetItem(item))
        self.ui.normalizeTable.setItem(idx, 4, QtWidgets.QTableWidgetItem(str(d.configuration.low_res_width)))
        item = QtWidgets.QTableWidgetItem(str(d.configuration.bck_position))
        item.setBackground(QtGui.QColor(200, 200, 200))
        self.ui.normalizeTable.setItem(idx, 5, QtWidgets.QTableWidgetItem(item))
        self.ui.normalizeTable.setItem(idx, 6, QtWidgets.QTableWidgetItem(str(d.configuration.bck_width)))
        self.main_window.auto_change_active = False

    def active_data_changed(self):
        """
            Actions to be taken once the active data set has changed
        """
        # If we update an entry, it's because that data is currently active.
        # Highlight it and un-highlight the other ones.
        self.main_window.auto_change_active = True
        idx = self._data_manager.find_active_data_id()
        for i in range(self.ui.reductionTable.rowCount()):
            item = self.ui.reductionTable.item(i, 0)
            if item is not None:
                if i == idx:
                    item.setBackground(QtGui.QColor(246, 213, 16))
                else:
                    item.setBackground(QtGui.QColor(255, 255, 255))

        idx = self._data_manager.find_active_direct_beam_id()
        for i in range(self.ui.normalizeTable.rowCount()):
            item = self.ui.normalizeTable.item(i, 0)
            if item is not None:
                if i == idx:
                    item.setBackground(QtGui.QColor(246, 213, 16))
                else:
                    item.setBackground(QtGui.QColor(255, 255, 255))
        self.main_window.auto_change_active = False

    def compute_offspec_on_change(self, force=False):
        """
            Compute off-specular as needed
        """
        prog = self.new_progress_reporter()
        has_changed_values = self.check_region_values_changed()
        offspec_data_exists = self._data_manager.is_offspec_available()
        logging.info("Exists %s %s", has_changed_values, offspec_data_exists)
        if force or has_changed_values>=0 or not offspec_data_exists:
            logging.info("Updating....")
            config = self.get_configuration()
            self._data_manager.update_configuration(configuration=config, active_only=False)
            self._data_manager.reduce_offspec(progress=prog)

    def compute_gisans_on_change(self, force=False, active_only=True):
        """
            Compute GISANS as needed
        """
        prog = self.new_progress_reporter()
        has_changed_values = self.check_region_values_changed()
        gisans_data_exists = self._data_manager.is_gisans_available(active_only=active_only)
        logging.info("Exists %s %s %s", force, has_changed_values, gisans_data_exists)
        if force or has_changed_values>=0 or not gisans_data_exists:
            logging.info("Updating....")
            config = self.get_configuration()
            self._data_manager.update_configuration(configuration=config, active_only=False)
            if active_only:
                self._data_manager.calculate_gisans(progress=prog)
            else:
                self._data_manager.reduce_gisans(active_only=active_only, progress=prog)

    def check_region_values_changed(self):
        """
            Return true if any of the parameters tied to a particular slot
            has changed.

            Some parameters are tied to the changeRegionValues() slot.
            There are time-consuming actions that we only want to take
            if those values actually changed, as opposed to the use simply
            clicking outside the box.

            Some parameters don't require a recalculation but simply a
            refreshing of the plots. Those are parameters such as scaling
            factors or the number of points clipped.

            Return values:
                -1 = no valid change
                 0 = replot needed
                 1 = recalculation needed
        """
        if self._data_manager.active_channel is None:
            return -1

        configuration = self._data_manager.active_channel.configuration
        valid_change = False
        replot_change = False

        # ROI parameters
        x_pos = self.ui.refXPos.value()
        x_width = self.ui.refXWidth.value()
        y_pos = self.ui.refYPos.value()
        y_width = self.ui.refYWidth.value()
        bck_pos = self.ui.bgCenter.value()
        bck_width = self.ui.bgWidth.value()

        valid_change = valid_change or \
            not configuration.peak_position == x_pos or \
            not configuration.peak_width == x_width

        valid_change = valid_change or \
            not configuration.low_res_position == y_pos or \
            not configuration.low_res_width == y_width

        valid_change = valid_change or \
            not configuration.bck_position == bck_pos or \
            not configuration.bck_width == bck_width

        try:
            scale = math.pow(10.0, self.ui.refScale.value())
        except:
            scale = 1
        replot_change = replot_change or \
            not configuration.scaling_factor == scale

        replot_change = replot_change or \
            not configuration.cut_first_n_points == self.ui.rangeStart.value()

        replot_change = replot_change or \
            not configuration.cut_last_n_points == self.ui.rangeEnd.value()

        valid_change = valid_change or \
            not configuration.subtract_background == self.ui.bgActive.isChecked()

        valid_change = valid_change or \
            not configuration.use_constant_q == self.ui.fanReflectivity.isChecked()

        valid_change = valid_change or \
            not configuration.use_dangle == self.ui.trustDANGLE.isChecked()

        valid_change = valid_change or \
            not configuration.set_direct_pixel == self.ui.set_dirpix_checkbox.isChecked()

        valid_change = valid_change or \
            not configuration.set_direct_angle_offset == self.ui.set_dangle0_checkbox.isChecked()

        if configuration.set_direct_pixel:
            valid_change = valid_change or \
                not configuration.direct_pixel_overwrite == self.ui.directPixelOverwrite.value()

        if configuration.set_direct_angle_offset:
            valid_change = valid_change or \
                not configuration.direct_angle_offset_overwrite == self.ui.dangle0Overwrite.value()

        # Final rebin
        valid_change = valid_change or \
            not configuration.do_final_rebin == self.ui.final_rebin_checkbox.isChecked()

        valid_change = valid_change or \
            not configuration.final_rebin_step == self.ui.q_rebin_spinbox.value()

        if valid_change:
            return 1
        if replot_change:
            return 0
        return -1

    def get_configuration(self):
        """
            Gather the reduction options.
        """
        if self._data_manager.active_channel is not None:
            configuration = self._data_manager.active_channel.configuration
        else:
            configuration = Configuration(self.main_window.settings)
        configuration.tof_bins = self.ui.eventTofBins.value()
        configuration.tof_bin_type = self.ui.eventBinMode.currentIndex()
        configuration.use_roi = self.ui.use_roi_checkbox.isChecked()
        configuration.update_peak_range = self.ui.fit_within_roi_checkbox.isChecked()
        configuration.use_roi_bck = self.ui.use_bck_roi_checkbox.isChecked()

        # Default ranges, using the current values
        x_pos = self.ui.refXPos.value()
        x_width = self.ui.refXWidth.value()
        y_pos = self.ui.refYPos.value()
        y_width = self.ui.refYWidth.value()
        bck_pos = self.ui.bgCenter.value()
        bck_width = self.ui.bgWidth.value()

        configuration.peak_position = x_pos
        configuration.peak_width = x_width
        configuration.low_res_position = y_pos
        configuration.low_res_width = y_width
        configuration.bck_position = bck_pos
        configuration.bck_width = bck_width

        configuration.force_peak_roi = not self.ui.actionAutomaticXPeak.isChecked()
        configuration.force_low_res_roi = not self.ui.actionAutoYLimits.isChecked()
        configuration.match_direct_beam = self.ui.actionAutoNorm.isChecked()

        # Use background on each side of the peak
        configuration.use_tight_bck = self.ui.use_side_bck_checkbox.isChecked()
        configuration.bck_offset = self.ui.side_bck_width.value()

        # Other reduction options
        configuration.subtract_background = self.ui.bgActive.isChecked()
        try:
            scale = math.pow(10.0, self.ui.refScale.value())
        except:
            scale = 1
        configuration.scaling_factor = scale
        configuration.cut_first_n_points = self.ui.rangeStart.value()
        configuration.cut_last_n_points = self.ui.rangeEnd.value()
        configuration.normalize_to_unity = self.ui.normalize_to_unity_checkbox.isChecked()
        configuration.total_reflectivity_q_cutoff = self.ui.normalization_q_cutoff_spinbox.value()
        configuration.wl_bandwidth = self.ui.bandwidth_spinbox.value()

        configuration.use_constant_q = self.ui.fanReflectivity.isChecked()
        configuration.use_dangle = self.ui.trustDANGLE.isChecked()
        configuration.set_direct_pixel = self.ui.set_dirpix_checkbox.isChecked()
        configuration.set_direct_angle_offset = self.ui.set_dangle0_checkbox.isChecked()
        configuration.direct_pixel_overwrite = self.ui.directPixelOverwrite.value()
        configuration.direct_angle_offset_overwrite = self.ui.dangle0Overwrite.value()
        configuration.sample_size = self.ui.sample_size_spinbox.value()
        configuration.do_final_rebin = self.ui.final_rebin_checkbox.isChecked()
        configuration.final_rebin_step = self.ui.q_rebin_spinbox.value()

        # UI elements
        configuration.normalize_x_tof = self.ui.normalizeXTof.isChecked()
        configuration.x_wl_map = self.ui.xLamda.isChecked()
        configuration.angle_map = self.ui.tthPhi.isChecked()
        configuration.log_1d = self.ui.logarithmic_y.isChecked()
        configuration.log_2d = self.ui.logarithmic_colorscale.isChecked()

        # Off-specular options
        if self.ui.kizmkfzVSqz.isChecked():
            configuration.off_spec_x_axis = Configuration.DELTA_KZ_VS_QZ
        elif self.ui.qxVSqz.isChecked():
            configuration.off_spec_x_axis = Configuration.QX_VS_QZ
        else:
            configuration.off_spec_x_axis = Configuration.KZI_VS_KZF
        configuration.off_spec_slice = self.ui.offspec_slice_checkbox.isChecked()
        configuration.off_spec_slice_qz_min = self.ui.slice_qz_min_spinbox.value()
        configuration.off_spec_slice_qz_max = self.ui.slice_qz_max_spinbox.value()
        #try:
        #    qz_list = self.ui.offspec_qz_list_edit.text()
        #    if len(qz_list) > 0:
        #        configuration.off_spec_qz_list = [float(x) for x in self.ui.offspec_qz_list_edit.text().split(',')]
        #except:
        #    logging.error("Could not parse off_spec_qz_list: %s", configuration.off_spec_qz_list)
        configuration.off_spec_err_weight = self.ui.offspec_err_weight_checkbox.isChecked()
        configuration.off_spec_nxbins = self.ui.offspec_rebin_x_bins_spinbox.value()
        configuration.off_spec_nybins = self.ui.offspec_rebin_y_bins_spinbox.value()
        configuration.off_spec_x_min = self.ui.offspec_x_min_spinbox.value()
        configuration.off_spec_x_max = self.ui.offspec_x_max_spinbox.value()
        configuration.off_spec_y_min = self.ui.offspec_y_min_spinbox.value()
        configuration.off_spec_y_max = self.ui.offspec_y_max_spinbox.value()

        # Off-spec smoothing options
        configuration.apply_smoothing = self.ui.offspec_smooth_checkbox.isChecked()

        # GISANS options
        configuration.gisans_wl_min = self.ui.gisans_wl_min_spinbox.value()
        configuration.gisans_wl_max = self.ui.gisans_wl_max_spinbox.value()
        configuration.gisans_wl_npts = self.ui.gisans_wl_npts_spinbox.value()
        configuration.gisans_qz_npts = self.ui.gisans_qz_npts_spinbox.value()
        configuration.gisans_qy_npts = self.ui.gisans_qy_npts_spinbox.value()
        configuration.gisans_use_pf = self.ui.gisans_pf_radio.isChecked()

        # Make the changes persistent
        configuration.to_q_settings(self.main_window.settings)
        return configuration

    def populate_from_configuration(self, configuration=None):
        """
            Set reduction options in UI, usually after loading
            a reduced data set.
        """
        if configuration is None:
            configuration = Configuration(self.main_window.settings)

        self.ui.eventTofBins.setValue(configuration.tof_bins)
        self.ui.eventBinMode.setCurrentIndex(configuration.tof_bin_type)
        self.ui.use_roi_checkbox.setChecked(configuration.use_roi)
        self.ui.fit_within_roi_checkbox.setChecked(configuration.update_peak_range)
        self.ui.use_bck_roi_checkbox.setChecked(configuration.use_roi_bck)

        self.ui.actionAutomaticXPeak.setChecked(not configuration.force_peak_roi)
        self.ui.actionAutoYLimits.setChecked(not configuration.force_low_res_roi)

        # Update reduction parameters
        self.ui.refXPos.setValue(configuration.peak_position)
        self.ui.refXWidth.setValue(configuration.peak_width)

        self.ui.refYPos.setValue(configuration.low_res_position)
        self.ui.refYWidth.setValue(configuration.low_res_width)

        self.ui.bgCenter.setValue(configuration.bck_position)
        self.ui.bgWidth.setValue(configuration.bck_width)

        # Use background on each side of the peak
        self.ui.use_side_bck_checkbox.setChecked(configuration.use_tight_bck)
        self.ui.side_bck_width.setValue(configuration.bck_offset)

        # Subtract background
        self.ui.bgActive.setChecked(configuration.subtract_background)
        # Scaling factor
        try:
            scale = math.log10(configuration.scaling_factor)
        except:
            scale = 0.0
        self.ui.refScale.setValue(scale)
        # Cut first and last points
        self.ui.rangeStart.setValue(configuration.cut_first_n_points)
        self.ui.rangeEnd.setValue(configuration.cut_last_n_points)
        self.ui.normalize_to_unity_checkbox.setChecked(configuration.normalize_to_unity)
        self.ui.normalization_q_cutoff_spinbox.setValue(configuration.total_reflectivity_q_cutoff)
        self.ui.bandwidth_spinbox.setValue(configuration.wl_bandwidth)

        self.ui.fanReflectivity.setChecked(configuration.use_constant_q)
        self.ui.trustDANGLE.setChecked(configuration.use_dangle)
        self.ui.set_dirpix_checkbox.setChecked(configuration.set_direct_pixel)
        self.ui.set_dangle0_checkbox.setChecked(configuration.set_direct_angle_offset)
        self.ui.directPixelOverwrite.setValue(configuration.direct_pixel_overwrite)
        self.ui.dangle0Overwrite.setValue(configuration.direct_angle_offset_overwrite)
        self.ui.sample_size_spinbox.setValue(configuration.sample_size)
        self.ui.final_rebin_checkbox.setChecked(configuration.do_final_rebin)
        self.ui.q_rebin_spinbox.setValue(configuration.final_rebin_step)

        # UI elements
        self.ui.normalizeXTof.setChecked(configuration.normalize_x_tof)
        self.ui.xLamda.setChecked(configuration.x_wl_map)
        self.ui.tthPhi.setChecked(configuration.angle_map)
        self.ui.logarithmic_y.setChecked(configuration.log_1d)
        self.ui.logarithmic_colorscale.setChecked(configuration.log_2d)

        # Off-specular options
        if configuration.off_spec_x_axis == Configuration.DELTA_KZ_VS_QZ:
            self.ui.kizmkfzVSqz.setChecked(True)
        elif configuration.off_spec_x_axis == Configuration.QX_VS_QZ:
            self.ui.qxVSqz.setChecked(True)
        else:
            self.ui.kizVSkfz.setChecked(True)
        self.ui.offspec_slice_checkbox.setChecked(configuration.off_spec_slice)
        #self.ui.offspec_qz_list_edit.setText(','.join([str(x) for x in configuration.off_spec_qz_list]))
        self.ui.slice_qz_min_spinbox.setValue(configuration.off_spec_slice_qz_min)
        self.ui.slice_qz_max_spinbox.setValue(configuration.off_spec_slice_qz_max)
        self.ui.offspec_err_weight_checkbox.setChecked(configuration.off_spec_err_weight)
        self.ui.offspec_rebin_x_bins_spinbox.setValue(configuration.off_spec_nxbins)
        self.ui.offspec_rebin_y_bins_spinbox.setValue(configuration.off_spec_nybins)
        self.ui.offspec_x_min_spinbox.setValue(configuration.off_spec_x_min)
        self.ui.offspec_x_max_spinbox.setValue(configuration.off_spec_x_max)
        self.ui.offspec_y_min_spinbox.setValue(configuration.off_spec_y_min)
        self.ui.offspec_y_max_spinbox.setValue(configuration.off_spec_y_max)

        # Off-spec smoothing options
        self.ui.offspec_smooth_checkbox.setChecked(configuration.apply_smoothing)

        # GISANS options
        self.ui.gisans_wl_min_spinbox.setValue(configuration.gisans_wl_min)
        self.ui.gisans_wl_max_spinbox.setValue(configuration.gisans_wl_max)
        self.ui.gisans_wl_npts_spinbox.setValue(configuration.gisans_wl_npts)
        self.ui.gisans_qz_npts_spinbox.setValue(configuration.gisans_qz_npts)
        self.ui.gisans_qy_npts_spinbox.setValue(configuration.gisans_qy_npts)
        self.ui.gisans_pf_radio.setChecked(configuration.gisans_use_pf)

    def stitch_reflectivity(self):
        """
            Stitch the reflectivity parts and normalize to 1.
        """
        # Update the configuration so we can remember the cutoff value
        # later if it was changed
        self.get_configuration()
        self._data_manager.stitch_data_sets(normalize_to_unity=self.ui.normalize_to_unity_checkbox.isChecked(),
                                            q_cutoff=self.ui.normalization_q_cutoff_spinbox.value())

        for i in range(len(self._data_manager.reduction_list)):
            xs = self._data_manager.active_channel.name
            d = self._data_manager.reduction_list[i].cross_sections[xs]
            self.ui.reductionTable.setItem(i, 1,
                                           QtWidgets.QTableWidgetItem("%.4f"%(d.configuration.scaling_factor)))

        self.main_window.initiate_reflectivity_plot.emit(False)

    def trim_data_to_normalization(self):
        """
            Cut the start and end of the active data set to 5% of its
            maximum intensity.
        """
        trim_points = self._data_manager.get_trim_values()
        if trim_points is not None:
            self.ui.rangeStart.setValue(trim_points[0])
            self.ui.rangeEnd.setValue(trim_points[1])
            self.update_tables()
            self.main_window.initiate_reflectivity_plot.emit(False)
        else:
            self.report_message("No direct beam found to trim data", pop_up=False)

    def strip_overlap(self):
        """
            Remove overlapping points in the reflecitviy, cutting always from the lower Qz
            measurements.
        """
        self._data_manager.strip_overlap()

        for i in range(len(self._data_manager.reduction_list)):
            xs = self._data_manager.active_channel.name
            d = self._data_manager.reduction_list[i].cross_sections[xs]
            self.ui.reductionTable.setItem(i, 3,
                                           QtWidgets.QTableWidgetItem(str(d.configuration.cut_last_n_points)))

        self.main_window.initiate_reflectivity_plot.emit(False)

    def report_message(self, message, informative_message=None,
                       detailed_message=None, pop_up=False, is_error=False):
        """
            Report an error.
            :param str message: message string to be reported
            :param str informative_message: extra information
            :param str detailed_message: detailed message for the log
            :param bool pop_up: if True, a dialog will pop up
            :param bool is_error: if True, the message is logged on the error channel
        """
        self.status_message.setText(message)
        if is_error:
            logging.error(message)
            if detailed_message is not None:
                logging.error(detailed_message)
        elif pop_up:
            logging.warning(message)
        else:
            logging.info(message)

        if pop_up:
            msg = QtWidgets.QMessageBox()
            msg.setIcon(QtWidgets.QMessageBox.Warning)

            msg.setText(message)
            msg.setWindowTitle("Information")
            if informative_message is not None:
                msg.setInformativeText(informative_message)
            if detailed_message is not None:
                msg.setDetailedText(detailed_message)
            msg.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msg.exec_()

    def show_results(self):
        """
            Pop up the result viewer
        """
        from ..result_viewer import ResultViewer
        dialog=ResultViewer(self.main_window, self._data_manager)
        dialog.specular_compare_widget.ui.refl_preview_checkbox.setChecked(True)
        self.main_window.update_specular_viewer.connect(dialog.update_specular)
        self.main_window.update_off_specular_viewer.connect(dialog.update_off_specular)
        self.main_window.update_gisans_viewer.connect(dialog.update_gisans)
        dialog.show()

# -*- coding: utf-8 -*-
#pylint: disable=invalid-name, line-too-long, too-many-public-methods, too-many-instance-attributes, wrong-import-order, bare-except
"""
    Main application window
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
from PyQt5 import QtCore, QtWidgets
import reflectivity_ui
import reflectivity_ui.interfaces.generated.ui_main_window
from reflectivity_ui.interfaces.event_handlers.plot_handler import PlotHandler
from reflectivity_ui.interfaces.event_handlers.main_handler import MainHandler
from .data_manager import DataManager
from .plotting import PlotManager
from .reduction_dialog import ReductionDialog
from .event_handlers.progress_reporter import ProgressReporter

class MainWindow(QtWidgets.QMainWindow,
                 reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow):
    """
        Main application window
    """
    # UI events
    file_loaded_signal = QtCore.pyqtSignal()
    initiate_projection_plot = QtCore.pyqtSignal(bool)
    initiate_reflectivity_plot = QtCore.pyqtSignal(bool)
    _gisansThread = None

    def __init__(self):
        """
            Initialization
        """
        # Base class
        QtWidgets.QMainWindow.__init__(self)

        # Initialize the UI widgets
        self.ui = reflectivity_ui.interfaces.generated.ui_main_window.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowTitle(u'QuickNXS %s' % reflectivity_ui.__version__)

        # Application settings
        self.settings = QtCore.QSettings('.refredm')

        # Object managers
        self.data_manager = DataManager(self.settings.value('current_directory', os.path.expanduser('~')))
        self.plot_manager = PlotManager(self)

        self.auto_change_active=False

        # Event handlers
        self.plot_handler = PlotHandler(self)
        self.file_handler = MainHandler(self)
        self.ui.compare_widget.data_manager = self.data_manager

        # Initialization for specific instrument
        # Retrieve configuration from config and enable/disable features
        self.initialize_instrument()
        self.hide_unsupported()

        # UI events
        self.file_loaded_signal.connect(self.file_handler.update_info)
        self.file_loaded_signal.connect(self.file_handler.update_daslog)
        self.file_loaded_signal.connect(self.plotActiveTab)
        self.initiate_projection_plot.connect(self.plot_manager.plot_projections)

        self.initiate_reflectivity_plot.connect(self.plot_manager.plot_refl)

    def closeEvent(self, event):
        """ Close UI event """
        self.file_handler.get_configuration()
        event.accept()

    def keyPressEvent(self, event):
        """ UI event """
        if event.modifiers()==QtCore.Qt.ControlModifier:
            self.plot_handler.control_down=True
        else:
            self.plot_handler.control_down=False

    def keyReleaseEvent(self, event):
        """ UI event """
        self.plot_handler.control_down=False

    def initialize_instrument(self):
        """
            Initialize instrument according to the instrument
            and saved parameters
        """
        for i in range(1, 12):
            getattr(self.ui, 'selectedChannel%i'%i).hide()
        self.ui.selectedChannel0.show()
        self.ui.selectedChannel0.setText(u"None")

        self.file_handler.populate_from_configuration()

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
        self.ui.counts_roi_label.hide()
        self.ui.eventActive.hide()

    # Actions defined in Qt Designer
    def file_open_dialog(self):
        """
            Show a dialog to open a new file.
            TODO: consider multiple selection. In this case QuickNXS tries to automatically sort and reduce.
        """
        self.file_handler.file_open_dialog()

    def file_loaded(self):
        """
            Update UI after a file is loaded
        """
        self.file_handler.file_loaded()

    def file_open_from_list(self):
        """
            Called when a new file is selected from the file list.
        """
        if self.auto_change_active:
            return
        QtWidgets.QApplication.instance().processEvents()
        item=self.ui.file_list.currentItem()
        name=unicode(item.text())
        QtWidgets.QApplication.instance().processEvents()
        self.file_handler.open_file(os.path.join(self.data_manager.current_directory, name))

    def reload_file(self):
        """
            Reload the file that is currently selected form the list.
        """
        self.file_handler.open_file(self.data_manager.current_file, force=True)

    def changeActiveChannel(self):
        """
            The overview and reflectivity channel was changed. This
            recalculates already extracted reflectivities.
        """
        return self.file_loaded()

    def getNorm(self):
        """
            TODO: deal with this
            This is supposed to retrieve the normalization data for the active reflectivity
            data so that we can normalize the distributions we are plotting.
            See plotting.plot_xtof and plotting.plot_overview
        """
        return self.data_manager.get_active_direct_beam()

    def plotActiveTab(self):
        '''
            Select the appropriate function to plot all visible images.
        '''
        if self.data_manager.active_channel is None:
            return
        color = str(self.ui.color_selector.currentText())
        if color!=self.plot_manager.color and self.plot_manager.color is not None:
            self.plot_manager.color = color
            plots=[self.ui.xy_pp, self.ui.xy_mp, self.ui.xy_pm, self.ui.xy_mm,
                   self.ui.xtof_pp, self.ui.xtof_mp, self.ui.xtof_pm, self.ui.xtof_mm,
                   self.ui.xy_overview, self.ui.xtof_overview]
            for plot in plots:
                plot.clear_fig()
        elif self.plot_manager.color is None:
            self.plot_manager.color = color
        if self.ui.plotTab.currentIndex()!=4 and self._gisansThread:
            self._gisansThread.finished.disconnect()
            self._gisansThread.terminate()
            self._gisansThread.wait(100)
            self._gisansThread=None
        if self.ui.plotTab.currentIndex()==0:
            self.plot_manager.plot_overview()
        elif self.ui.plotTab.currentIndex()==1:
            self.plot_manager.plot_xy()
        elif self.ui.plotTab.currentIndex()==2:
            self.plot_manager.plot_xtof()
        elif self.ui.plotTab.currentIndex()==3:
            self.file_handler.compute_offspec_on_change()
            self.plot_manager.plot_offspec()
        elif self.ui.plotTab.currentIndex()==4:
            self.file_handler.compute_gisans_on_change(active_only=True)
            self.plot_manager.plot_gisans()
        elif self.ui.plotTab.currentIndex()==6:
            self.ui.compare_widget.draw()

    def toggleColorbars(self):
        """ Refresh plots because of a color or scale change """
        plots=[self.ui.xy_pp, self.ui.xy_mp, self.ui.xy_pm, self.ui.xy_mm,
               self.ui.xtof_pp, self.ui.xtof_mp, self.ui.xtof_pm, self.ui.xtof_mm,
               self.ui.xy_overview, self.ui.xtof_overview,
               self.ui.offspec_pp, self.ui.offspec_mp, self.ui.offspec_pm, self.ui.offspec_mm]
        for plot in plots:
            plot.clear_fig()
        self.plotActiveTab()

    def changeRegionValues(self):
        """
            Called when the reflectivity extraction region has been changed.
            Sets up a trigger to replot the reflectivity with a delay so
            a subsequent change can occur without several replots.
        """
        if self.auto_change_active:
            return

        change_type = self.file_handler.check_region_values_changed()
        if change_type >= 0:
            configuration = self.file_handler.get_configuration()

            if self.data_manager.active_channel is not None:
                active_only = not self.ui.action_use_common_ranges.isChecked()
                self.data_manager.update_configuration(configuration=configuration, active_only=active_only)
                self.plot_handler.change_region_values()
                self.file_handler.update_calculated_data()

                # Update the reduction/direct beam tables if this data set is in it
                self.file_handler.update_tables()

                QtWidgets.QApplication.instance().processEvents()
                if change_type > 0:
                    try:
                        self.data_manager.calculate_reflectivity(configuration=configuration, active_only=active_only)
                    except:
                        self.file_handler.report_message("There was a problem updating the reflectivity",
                                                         pop_up=False)
                        logging.error("There was a problem updating the reflectivity\n%s", sys.exc_value)
                self.plot_manager.plot_refl()

    def reductionTableChanged(self, item):
        '''
        Perform action upon change in data reduction list.
        '''
        self.file_handler.reduction_table_changed(item)

    def reduction_cell_activated(self, row, col):
        """
            Select a data set when the user double-clicks on a run number (col 0).
            in the reduction table.
            :param int row: row index
            :param int col: column index
        """
        if col == 0:
            self.data_manager.set_active_data_from_reduction_list(row)
            self.file_loaded()
            self.file_handler.active_data_changed()

    def direct_beam_cell_activated(self, row, col):
        """
            Select a data set when the user double-clicks on a run number (col 0).
            in the direct beam table.
            :param int row: row index
            :param int col: column index
        """
        if col == 0:
            self.data_manager.set_active_data_from_direct_beam_list(row)
            self.file_loaded()
            #self.file_handler.active_data_changed()

    def replotProjections(self):
        """ Signal handling """
        self.initiate_projection_plot.emit(True)
        self.initiate_reflectivity_plot.emit(True)

    def addRefList(self):
        """ Signal handling """
        self.file_handler.add_reflectivity()

    def removeRefList(self):
        """ Signal handling """
        self.file_handler.remove_reflectivity()

    def clearRefList(self):
        """ Signal handling """
        self.file_handler.clear_reflectivity()

    def setNorm(self):
        """ Signal handling """
        self.file_handler.add_direct_beam()

    def remove_normalization(self):
        """ Signal handling """
        self.file_handler.remove_direct_beam()

    def clearNormList(self):
        """ Signal handling """
        self.file_handler.clear_direct_beams()

    def match_direct_beam_clicked(self):
        """
            Find the best direct beam run for the activate data set
            and compute the reflectivity as needed.
        """
        if self.data_manager.find_best_direct_beam():
            self.file_handler.update_tables()
            self.file_handler.update_calculated_data()
            QtWidgets.QApplication.instance().processEvents()
            try:
                self.data_manager.calculate_reflectivity()
            except:
                self.file_handler.report_message("There was a problem updating the reflectivity",
                                                 pop_up=False)
                logging.error("There was a problem updating the reflectivity\n%s", sys.exc_value)
            self.initiate_reflectivity_plot.emit(True)

    def openByNumber(self):
        """ Signal handling """
        self.file_handler.open_run_number()

    def refresh_offspec(self):
        """
            Refresh / recalculate the off-specular plots
        """
        self.file_handler.compute_offspec_on_change(force=True)
        self.plot_manager.plot_offspec()

    def change_offspec_colorscale(self):
        """
            Change the intensity limits for the color scale of the off-specular plots
        """
        self.plot_handler.change_offspec_colorscale()

    def cutPoints(self):
        """
            Cut the start and end of the active data set to 5% of its
            maximum intensity.
        """
        self.file_handler.trim_data_to_normalization()

    def stripOverlap(self):
        """
            Remove overlapping points in the reflecitviy, cutting always from the lower Qz
            measurements.
        """
        self.file_handler.strip_overlap()

    def normalizeTotalReflection(self):
        """
            Stitch the reflectivity parts and normalize to 1.
        """
        self.file_handler.stitch_reflectivity()

    def autoRef(self):
        self.file_handler.automated_file_selection()

    def reduceDatasets(self):
        '''
        Open a dialog to select reduction options for the current list of
        reduction items.
        '''
        if len(self.data_manager.reduction_list)==0:
            self.file_handler.report_message("The data to be reduced must be added to the reduction table",
                                             pop_up=True)
            return
        dialog=ReductionDialog(self)
        dialog.exec_()
        output_options = dialog.get_options()
        dialog.destroy()

        if output_options is not None:
            configuration = self.file_handler.get_configuration()
            output_options['off_spec_x_axis'] = configuration.off_spec_x_axis
            output_options['off_spec_slice'] = configuration.off_spec_slice
            output_options['off_spec_qz_list'] = configuration.off_spec_qz_list
            output_options['off_spec_err_weight'] = configuration.off_spec_err_weight
            output_options['off_spec_nxbins'] = configuration.off_spec_nxbins
            output_options['off_spec_nybins'] = configuration.off_spec_nybins

            from .data_handling.processing_workflow import ProcessingWorkflow
            wrk = ProcessingWorkflow(self.data_manager, output_options)
            wrk.execute(self.file_handler.new_progress_reporter())

    def loadExtraction(self):
        self.file_handler.open_reduced_file_dialog()

    def refresh_file_list(self):
        self.file_handler.update_file_list()

    # Un-used UI signals
    #pylint: disable=missing-docstring, multiple-statements, no-self-use
    def change_gisans_colorscale(self): return NotImplemented
    def fileOpenSumDialog(self): return NotImplemented

    # From the Advanced menu
    def open_advanced_background(self): return NotImplemented
    def open_polarization_window(self): return NotImplemented
    def open_rawdata_dialog(self): return NotImplemented

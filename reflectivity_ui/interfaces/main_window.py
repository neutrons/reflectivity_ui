# -*- coding: utf-8 -*-
#pylint: disable=invalid-name, line-too-long, too-many-public-methods, too-many-instance-attributes
"""
    Main application window
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
import copy
from PyQt5 import QtGui, QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_main_window

from reflectivity_ui.interfaces.event_handlers.plot_handler import PlotHandler
from reflectivity_ui.interfaces.event_handlers.main_handler import MainHandler
from .data_manager import DataManager
from .plotting import PlotManager

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

        # Application settings
        self.settings = QtCore.QSettings('.refredm')

        # Object managers
        self.data_manager = DataManager(self.settings.value('current_directory', os.path.expanduser('~')))
        self.plot_manager = PlotManager(self)

        self.auto_change_active=False

        # Event handlers
        self.plot_handler = PlotHandler(self)
        self.file_handler = MainHandler(self)

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
        self.file_handler.get_configuration()
        event.accept()

    def keyPressEvent(self, event):
        if event.modifiers()==QtCore.Qt.ControlModifier:
            self.plot_handler.control_down=True
        else:
            self.plot_handler.control_down=False

    def keyReleaseEvent(self, event):
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
        self.ui.eventActive.hide()

        # Hide quick reduce button
        self.ui.load_live_data_button.hide()

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
        return None

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
        if self.ui.plotTab.currentIndex()!=4 and self._gisansThread: #TODO
            self._gisansThread.finished.disconnect()
            self._gisansThread.terminate()
            self._gisansThread.wait(100)
            self._gisansThread=None
        if self.ui.plotTab.currentIndex()==0:
            self.plot_manager.plot_overview()
        if self.ui.plotTab.currentIndex()==1:
            self.plot_manager.plot_xy()
        if self.ui.plotTab.currentIndex()==2:
            self.plot_manager.plot_xtof()
        if self.ui.plotTab.currentIndex()==3:
            self.plot_offspec() #TODO
        if self.ui.plotTab.currentIndex()==4:
            self.plot_gisans() #TODO

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

        if self.file_handler.check_region_values_changed():
            configuration = self.file_handler.get_configuration()

            if self.data_manager.active_channel is not None:
                active_only = not self.ui.action_use_common_ranges.isChecked()
                self.data_manager.update_configuration(configuration=configuration, active_only=active_only)
                self.plot_handler.change_region_values()
                self.file_handler.update_calculated_data()

                # Update the reduction/direct beam tables if this data set is in it
                self.file_handler.update_tables()

                QtWidgets.QApplication.instance().processEvents()
                try:
                    self.data_manager.calculate_reflectivity(configuration=configuration, active_only=active_only)
                except:
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
            self.file_handler.active_data_changed()

    def replotProjections(self):
        self.initiate_projection_plot.emit(True)
        self.initiate_reflectivity_plot.emit(True)

    def addRefList(self):
        self.file_handler.add_reflectivity()
    def removeRefList(self):
        self.file_handler.remove_reflectivity()
    def clearRefList(self):
        self.file_handler.clear_reflectivity()

    def setNorm(self):
        self.file_handler.add_direct_beam()

    def remove_normalization(self):
        self.file_handler.remove_direct_beam()

    def clearNormList(self):
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
            self.data_manager.calculate_reflectivity()
            self.initiate_reflectivity_plot.emit(True)

    def openByNumber(self):
        self.file_handler.open_run_number()

    def normalizeTotalReflection(self): return NotImplemented
    def reduceDatasets(self): return NotImplemented
    def loadExtraction(self): return NotImplemented
    def change_offspec_colorscale(self): return NotImplemented
    def change_gisans_colorscale(self): return NotImplemented
    def clip_offspec_colorscale(self): return NotImplemented
    def fileOpenSumDialog(self): return NotImplemented
    def overwriteChanged(self): return NotImplemented
    def cutPoints(self): return NotImplemented
    def autoRef(self): return NotImplemented
    def stripOverlap(self): return NotImplemented
    def live_open(self): return NotImplemented

    # From the Advanced menu
    def overwriteDirectBeam(self): return NotImplemented
    def open_advanced_background(self): return NotImplemented
    def clearOverwrite(self): return NotImplemented
    def open_polarization_window(self): return NotImplemented
    def open_rawdata_dialog(self): return NotImplemented

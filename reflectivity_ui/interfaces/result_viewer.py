"""
   Dialog to show final reduced data.
"""
#pylint: disable=bare-except

import logging
import numpy as np
from PyQt5 import QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_result_viewer
import reflectivity_ui.interfaces.generated.mplwidget as mpl

class ResultViewer(QtWidgets.QDialog, reflectivity_ui.interfaces.generated.ui_result_viewer.Ui_Dialog):
    """
        Reduction dialog
    """
    default_template = '(instrument)_{numbers}_{item}_{state}.{type}'

    def __init__(self, parent, data_manager):
        super(ResultViewer, self).__init__(parent)
        self.setupUi(self)
        self.resize(1024, 1024)
        self.settings = QtCore.QSettings('.refredm')
        self.data_manager = data_manager
        self.main_window = parent
        self.specular_compare_widget.data_manager = self.data_manager
        self._gisans_reference = None

    def update_active_tab(self):
        if self.tabWidget.currentIndex()==0:
            self.update_specular()
        elif self.tabWidget.currentIndex()==1:
            self.update_off_specular()
        elif self.tabWidget.currentIndex()==2:
            self.update_gisans()

    def update_specular(self):
        self.specular_compare_widget.refl_preview()

    def update_off_specular(self, crop=False):
        """
            Update the result viewer with the latest off-specular calculations.
            :param bool crop: if True, all the plots will be cropped to the ++ cross-section
        """
        off_spec_data = self.data_manager.cached_offspec
        if off_spec_data is None:
            return

        xlim = None
        ylim = None
        if crop and self.offspec_pp_plot.cplot is not None:
            xlim = self.offspec_pp_plot.canvas.ax.get_xlim()
            ylim = self.offspec_pp_plot.canvas.ax.get_ylim()

        data_set_keys = list(self.data_manager.data_sets.keys())

        if len(data_set_keys) > 4:
            logging.error("Too many cross-sections for plotting: %s", str(len(data_set_keys)))

        plots=[self.offspec_pp_plot, self.offspec_mm_plot,
               self.offspec_pm_plot, self.offspec_mp_plot]
        for plot in plots:
            plot.clear()

        for i in range(len(data_set_keys), 4):
            if plots[i].cplot is not None:
                plots[i].draw()
            plots[i].hide()

        i_min=10**self.offspec_intensity_min.value()
        i_max=10**self.offspec_intensity_max.value()

        for i, channel in enumerate(data_set_keys):
            plot = plots[i]
            plot.show()
            plots[i].clear_fig()
            _data = off_spec_data[channel][0].T
            plots[i].pcolormesh(_data[0], _data[1], _data[2], log=True,
                                imin=i_min, imax=i_max)
            plots[i].set_xlabel('%s [%s]' % (off_spec_data['columns'][0], off_spec_data['units'][0]))
            plots[i].set_ylabel('%s [%s]' % (off_spec_data['columns'][1], off_spec_data['units'][1]))
            plots[i].set_title(channel)
            if plots[i].cplot is not None:
                plots[i].cplot.set_clim([i_min, i_max])
                if xlim is not None and ylim is not None:
                    plots[i].canvas.ax.set_xlim(*xlim)
                    plots[i].canvas.ax.set_ylim(*ylim)
            plots[i].draw()

    def apply_offspec_crop(self):
        self.update_off_specular(crop=True)

    def reset_offspec_crop(self):
        self.update_off_specular(crop=False)

    def apply_gisans_crop(self):
        self.update_gisans(crop=True)

    def reset_gisans_crop(self):
        self.update_gisans(crop=False)

    def _plot_gisans(self, gisans_data, channel, layout, i_min, i_max, xlim=None, ylim=None):
        _data = gisans_data[channel][0].T
        gisans_plot = mpl.MPLWidget(self)
        gisans_plot.setMinimumSize(QtCore.QSize(0, 250))
        if self._gisans_reference is None:
            self._gisans_reference = gisans_plot
        layout.addWidget(gisans_plot)#, i_row, i_col)
        
        gisans_plot.pcolormesh(_data[0], _data[1], _data[2], log=True,
                            imin=i_min, imax=i_max)
        gisans_plot.set_xlabel('%s [%s]' % (gisans_data['columns'][0], gisans_data['units'][0]))
        gisans_plot.set_ylabel('%s [%s]' % (gisans_data['columns'][1], gisans_data['units'][1]))
        gisans_plot.set_title(channel)
        if xlim is not None and ylim is not None:
            gisans_plot.canvas.ax.set_xlim(*xlim)
            gisans_plot.canvas.ax.set_ylim(*ylim)

        gisans_plot.draw()
        return gisans_plot

    def update_gisans(self, crop=False):
        """
            Update the results viewer with the latest GISANS calculations
            :param bool crop: if True, all the plots will be cropped to the ++ cross-section
        """
        logging.info("Updating GISANS")
        gisans_data = self.data_manager.cached_gisans
        if gisans_data is None:
            logging.info("Nothing to plot for GISANS")
            return

        # Clear everything
        xlim = None
        ylim = None
        if crop and self._gisans_reference is not None and self._gisans_reference.cplot is not None:
            xlim = self._gisans_reference.canvas.ax.get_xlim()
            ylim = self._gisans_reference.canvas.ax.get_ylim()

        self._gisans_reference = None

        data_set_keys = list(self.data_manager.data_sets.keys())
        if len(data_set_keys) > 4:
            logging.error("Too many cross-sections for plotting: %s", str(len(data_set_keys)))

        layouts = [self.gisans_pp_layout, self.gisans_mm_layout,
                   self.gisans_pm_layout, self.gisans_mp_layout]

        i_min=10**self.gisans_intensity_min.value()
        i_max=10**self.gisans_intensity_max.value()

        for i, pol_state in enumerate(data_set_keys):
            clear_layout(layouts[i])
            logging.info("State: %s" % pol_state)
            for j, channel in enumerate(gisans_data["cross_section_bins"][pol_state]):
                logging.info("  channel: %s", channel)
                _plot = self._plot_gisans(gisans_data, channel, layouts[i], i_min, i_max, xlim=xlim, ylim=ylim)
                if i == 0 and j == 0:
                    self._gisans_reference = _plot

def clear_layout(layout):
    while layout.count() > 0:
        item = layout.takeAt(0)
        if not item:
            continue
        w = item.widget()
        if w:
            w.deleteLater()

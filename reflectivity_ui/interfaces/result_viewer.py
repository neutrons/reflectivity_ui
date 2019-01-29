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

    def __init__(self, parent, data_manager):
        super(ResultViewer, self).__init__(parent)
        self.setupUi(self)
        self.resize(1024, 1024)
        self.settings = QtCore.QSettings('.refredm')
        self.data_manager = data_manager
        self.main_window = parent
        self.specular_compare_widget.data_manager = self.data_manager

    def update_active_tab(self):
        if self.tabWidget.currentIndex()==0:
            self.update_specular()
        elif self.tabWidget.currentIndex()==1:
            self.update_off_specular()
        elif self.tabWidget.currentIndex()==2:
            # Not yet implemented
            pass

    def update_specular(self):
        self.specular_compare_widget.refl_preview()

    def update_off_specular(self, crop=False):
        off_spec_data = self.data_manager.cached_offspec
        if off_spec_data is None:
            return

        xlim = None
        ylim = None
        if crop and self.offspec_pp_plot.cplot is not None:
            xlim = self.offspec_pp_plot.canvas.ax.get_xlim()
            ylim = self.offspec_pp_plot.canvas.ax.get_ylim()

        data_set_keys = off_spec_data['cross_sections'].keys()

        plots=[self.offspec_pp_plot, self.offspec_mm_plot,
               self.offspec_pm_plot, self.offspec_mp_plot]
        for plot in plots:
            plot.clear()

        for i in range(len(data_set_keys), 4):
            if plots[i].cplot is not None:
                plots[i].draw()

        if len(data_set_keys)>1:
            self.main_window.ui.offspec_mm.show()
            if len(data_set_keys)==4:
                self.offspec_mp_plot.show()
                self.offspec_pm_plot.show()
            else:
                self.offspec_mp_plot.hide()
                self.offspec_pm_plot.hide()
        else:
            self.offspec_mp_plot.hide()
            self.offspec_pm_plot.hide()
            self.offspec_mm_plot.hide()

        i_min=10**self.offspec_intensity_min.value()
        i_max=10**self.offspec_intensity_max.value()

        for i, channel in enumerate(data_set_keys):
            plot = plots[i]
            plots[i].clear_fig()
            _data = off_spec_data[channel][0].T
            plots[i].pcolormesh(_data[0], _data[1], _data[2], log=True,
                                imin=i_min, imax=i_max,
                                shading='gouraud')
            plots[i].set_xlabel(u'%s [%s]' % (off_spec_data['columns'][0], off_spec_data['units'][0]))
            plots[i].set_ylabel(u'%s [%s]' % (off_spec_data['columns'][1], off_spec_data['units'][1]))
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

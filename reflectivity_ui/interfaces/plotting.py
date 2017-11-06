# -*- coding: utf-8 -*-
import numpy as np
import logging

class PlotManager(object):
    def __init__(self, main_window):
        self.main_window = main_window
        self.overview_lines=None
        self._x_projection = None
        self._y_projection = None
        self.proj_lines = None
        self.y_bg=0.
        self.color = None
        
        # Selection lines
        self.xy_x1 = None
        self.xy_x2 = None
        self.xy_y1 = None
        self.xy_y2 = None

        self.xtof_x1 = None
        self.xtof_x2 = None
        self.xtof_bck1 = None
        self.xtof_bck2 = None

    def plot_overview(self):
        '''
        X vs. Y and X vs. Tof for main channel.
        '''
        self.xy_x1 = None
        self.xy_x2 = None
        self.xy_y1 = None
        self.xy_y2 = None
        self.xtof_x1 = None
        self.xtof_x2 = None
        self.xtof_bck1 = None
        self.xtof_bck2 = None
        main_window = self.main_window
        main_window.ui.xy_overview.clear()
        main_window.ui.xtof_overview.clear()

        data=main_window.data_manager.active_channel
        if data.total_counts==0:
            main_window.ui.xy_overview.draw()
            main_window.ui.xtof_overview.draw()
            return

        xy=data.xydata
        xtof=data.xtofdata
        ref_norm=main_window.getNorm()
        if main_window.ui.normalizeXTof.isChecked() and ref_norm is not None:
            ref_norm=ref_norm.Rraw
            # normalize ToF dataset for wavelength distribution
            ref_norm=np.where(ref_norm>0., ref_norm, 1.)
            xtof=xtof.astype(float)/ref_norm[np.newaxis, :]
        xy_imin=xy[xy>0].min()
        xy_imax=xy.max()
        tof_imin=xtof[xtof>0].min()
        tof_imax=xtof.max()

        # for lines of the current extraction area
        x_peak=main_window.ui.refXPos.value()
        x_width=main_window.ui.refXWidth.value()
        y_pos=main_window.ui.refYPos.value()
        y_width=main_window.ui.refYWidth.value()
        bg_pos=main_window.ui.bgCenter.value()
        bg_width=main_window.ui.bgWidth.value()

        # XY plot
        if main_window.ui.tthPhi.isChecked():
            rad_per_pixel=data.det_size_x/data.dist_sam_det/data.xydata.shape[1]
            phi_range=xy.shape[0]*rad_per_pixel*180./np.pi
            tth_range=xy.shape[1]*rad_per_pixel*180./np.pi
            phi0=main_window.ui.refYPos.value()*rad_per_pixel*180./np.pi
            tth0=(data.dangle-data.dangle0)-(xy.shape[1]-data.dpix)*rad_per_pixel*180./np.pi
            main_window.ui.xy_overview.imshow(xy, log=main_window.ui.logarithmic_colorscale.isChecked(),
                                              aspect='auto', cmap=self.color, origin='lower',
                                              extent=[tth_range+tth0, tth0, phi0, phi0-phi_range])
            main_window.ui.xy_overview.set_xlabel(u'2$\\Theta{}$ [°]')
            main_window.ui.xy_overview.set_ylabel(u'$\\phi{}$ [°]')
            main_window.ui.xy_overview.cplot.set_clim([xy_imin, xy_imax])
        else:
            main_window.ui.xy_overview.imshow(xy, log=main_window.ui.logarithmic_colorscale.isChecked(),
                                              aspect='auto', cmap=self.color, origin='lower')
            main_window.ui.xy_overview.set_xlabel(u'x [pix]')
            main_window.ui.xy_overview.set_ylabel(u'y [pix]')
            main_window.ui.xy_overview.cplot.set_clim([xy_imin, xy_imax])

            if self.xy_x1 is None:
                self.xy_x1 = main_window.ui.xy_overview.canvas.ax.axvline(x_peak-x_width/2., color='#aa0000')
                self.xy_x2 = main_window.ui.xy_overview.canvas.ax.axvline(x_peak+x_width/2., color='#aa0000')
                self.xy_y1 = main_window.ui.xy_overview.canvas.ax.axhline(y_pos-y_width/2., color='#00aa00')
                self.xy_y2 = main_window.ui.xy_overview.canvas.ax.axhline(y_pos+y_width/2., color='#00aa00')
            else:
                self.xy_x1.set_xdata([x_peak-x_width/2., x_peak-x_width/2.])
                self.xy_x2.set_xdata([x_peak+x_width/2., x_peak+x_width/2.])
                self.xy_y1.set_ydata([y_pos-y_width/2., y_pos-y_width/2.])
                self.xy_y2.set_ydata([y_pos+y_width/2., y_pos+y_width/2.])

        # XToF plot
        if main_window.ui.xLamda.isChecked():
            main_window.ui.xtof_overview.imshow(xtof[::-1], log=main_window.ui.logarithmic_colorscale.isChecked(),
                                         aspect='auto', cmap=self.color,
                                         extent=[data.wavelength[0], data.wavelength[-1], 0, data.x.shape[0]-1])
            main_window.ui.xtof_overview.set_xlabel(u'$\\lambda{}$ [Å]')
        else:
            main_window.ui.xtof_overview.imshow(xtof[::-1], log=main_window.ui.logarithmic_colorscale.isChecked(),
                                         aspect='auto', cmap=self.color,
                                         extent=[data.tof[0]*1e-3, data.tof[-1]*1e-3, 0, data.x.shape[0]-1])
            main_window.ui.xtof_overview.set_xlabel(u'ToF [ms]')
        main_window.ui.xtof_overview.set_ylabel(u'x [pix]')

        if self.xtof_x1 is None:
            self.xtof_x1 = main_window.ui.xtof_overview.canvas.ax.axhline(x_peak-x_width/2., color='#aa0000')
            self.xtof_x2 = main_window.ui.xtof_overview.canvas.ax.axhline(x_peak+x_width/2., color='#aa0000')
            self.xtof_bck1 = main_window.ui.xtof_overview.canvas.ax.axhline(bg_pos-bg_width/2., color='black')
            self.xtof_bck2 = main_window.ui.xtof_overview.canvas.ax.axhline(bg_pos+bg_width/2., color='black')
        else:
            self.xtof_x1.set_ydata([x_peak-x_width/2., x_peak-x_width/2.])
            self.xtof_x2.set_ydata([x_peak+x_width/2., x_peak+x_width/2.])
            self.xtof_bck1.set_ydata([bg_pos-bg_width/2., bg_pos-bg_width/2.])
            self.xtof_bck2.set_ydata([bg_pos+bg_width/2., bg_pos+bg_width/2.])

        main_window.ui.xtof_overview.cplot.set_clim([tof_imin, tof_imax])

        if main_window.ui.show_colorbars.isChecked() and main_window.ui.xy_overview.cbar is None:
            main_window.ui.xy_overview.cbar=main_window.ui.xy_overview.canvas.fig.colorbar(main_window.ui.xy_overview.cplot)
            main_window.ui.xtof_overview.cbar=main_window.ui.xtof_overview.canvas.fig.colorbar(main_window.ui.xtof_overview.cplot)
        main_window.ui.xy_overview.draw()
        main_window.ui.xtof_overview.draw()

    def plot_xy(self):
        """
            X vs. Y plots for all channels.
        """
        main_window = self.main_window
        data_set_keys = main_window.data_manager.data_sets.keys()
        plots=[main_window.ui.xy_pp, main_window.ui.xy_mm, main_window.ui.xy_pm, main_window.ui.xy_mp]
        for i in range(len(data_set_keys), 4):
            if plots[i].cplot is not None:
                plots[i].clear()
                plots[i].draw()
        imin=1e20
        imax=1e-20
        xynormed=[]

        for key in data_set_keys[:4]:
            dataset = main_window.data_manager.data_sets[key]
            d=dataset.xydata/dataset.proton_charge
            xynormed.append(d)
            if dataset.total_counts==0:
                continue
            imin=min(imin, d[d>0].min())
            imax=max(imax, d.max())

        if len(xynormed)>1:
            main_window.ui.frame_xy_mm.show()
            if len(xynormed)==4:
                main_window.ui.frame_xy_sf.show()
            else:
                main_window.ui.frame_xy_sf.hide()
        else:
            main_window.ui.frame_xy_mm.hide()
            main_window.ui.frame_xy_sf.hide()
    
        for i, datai in enumerate(xynormed):
            if main_window.data_manager.data_sets[data_set_keys[i]].total_counts==0:
                continue
            if main_window.ui.tthPhi.isChecked():
                plots[i].clear()
                rad_per_pixel=dataset.det_size_x/dataset.dist_sam_det/dataset.xydata.shape[1]
                phi_range=datai.shape[0]*rad_per_pixel*180./np.pi
                tth_range=datai.shape[1]*rad_per_pixel*180./np.pi
                phi0=main_window.ui.refYPos.value()*rad_per_pixel*180./np.pi
                tth0=(dataset.dangle-dataset.dangle0)-(datai.shape[1]-dataset.dpix)*rad_per_pixel*180./np.pi
                plots[i].imshow(datai, log=main_window.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, origin='lower',
                                extent=[tth_range+tth0, tth0, phi0, phi0-phi_range])
                plots[i].set_xlabel(u'2$\\Theta{}$ [°]')
                plots[i].set_ylabel(u'$\\phi{}$ [°]')
            else:
                plots[i].imshow(datai, log=main_window.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, origin='lower')
                plots[i].set_xlabel(u'x [pix]')
                plots[i].set_ylabel(u'y [pix]')
            plots[i].set_title(data_set_keys[i])
            if plots[i].cplot is not None:
                plots[i].cplot.set_clim([imin, imax])
            if plots[i].cplot is not None and main_window.ui.show_colorbars.isChecked() and plots[i].cbar is None:
                plots[i].cbar=plots[i].canvas.fig.colorbar(plots[i].cplot)
            plots[i].draw()

    def plot_xtof(self):
        """
            X vs. ToF plots for all channels.
        """
        main_window = self.main_window
        data_set_keys = main_window.data_manager.data_sets.keys()
        imin=1e20
        imax=1e-20
        xtofnormed=[]
        ref_norm=main_window.getNorm()
        if ref_norm is not None:
            ref_norm=ref_norm.Rraw
            ref_norm=np.where(ref_norm>0, ref_norm, 1.)

        for key in data_set_keys[:4]:
            dataset = main_window.data_manager.data_sets[key]
            d=dataset.xtofdata/dataset.proton_charge
            if main_window.ui.normalizeXTof.isChecked() and ref_norm is not None:
                # normalize all datasets for wavelength distribution
                d=d/ref_norm[np.newaxis, :]
            xtofnormed.append(d)
            if dataset.total_counts==0:
                continue
            imin=min(imin, d[d>0].min())
            imax=max(imax, d.max())
        wavelength = main_window.data_manager.active_channel.wavelength
        tof = main_window.data_manager.active_channel.tof
    
        plots=[main_window.ui.xtof_pp, main_window.ui.xtof_mm, main_window.ui.xtof_pm, main_window.ui.xtof_mp]
        for i in range(len(data_set_keys), 4):
            if plots[i].cplot is not None:
                plots[i].clear()
                plots[i].draw()

        if len(xtofnormed)>1:
            main_window.ui.frame_xtof_mm.show()
            if len(xtofnormed)==4:
                main_window.ui.frame_xtof_sf.show()
            else:
                main_window.ui.frame_xtof_sf.hide()
        else:
            main_window.ui.frame_xtof_mm.hide()
            main_window.ui.frame_xtof_sf.hide()
        for i, datai in enumerate(xtofnormed):
            if main_window.data_manager.data_sets[data_set_keys[i]].total_counts==0:
                continue
            if main_window.ui.xLamda.isChecked():
                plots[i].imshow(datai[::-1], log=main_window.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, extent=[wavelength[0], wavelength[-1], 0, datai.shape[0]-1])
                plots[i].set_xlabel(u'$\\lambda{}$ [Å]')
            else:
                plots[i].imshow(datai[::-1], log=main_window.ui.logarithmic_colorscale.isChecked(), imin=imin, imax=imax,
                                aspect='auto', cmap=self.color, extent=[tof[0]*1e-3, tof[-1]*1e-3, 0, datai.shape[0]-1])
                plots[i].set_xlabel(u'ToF [ms]')
            plots[i].set_title(data_set_keys[i])
            plots[i].set_ylabel(u'x [pix]')
            if plots[i].cplot is not None:
                plots[i].cplot.set_clim([imin, imax])
            if plots[i].cplot is not None and main_window.ui.show_colorbars.isChecked() and plots[i].cbar is None:
                plots[i].cbar=plots[i].canvas.fig.colorbar(plots[i].cplot)
            plots[i].draw()

    def plot_projections(self, preserve_lim=False):
        """
            Create projections of the data on the x and y axes.
            The x-projection can also be done be means of quantile calculation,
            which means that the ToF intensities are calculation which are
            exceeded by a certain number of points. This can be helpful to better
            separate the specular reflection from bragg-sheets
        """
        main_window = self.main_window
        if main_window.data_manager.active_channel is None:
            return
        data = main_window.data_manager.active_channel
        if data.total_counts == 0:
            main_window.ui.x_project.clear()
            main_window.ui.x_project.draw()
            main_window.ui.y_project.clear()
            main_window.ui.y_project.draw()
            return

        xproj=data.xdata
        yproj=data.ydata

        x_peak=main_window.ui.refXPos.value()
        x_width=main_window.ui.refXWidth.value()
        y_pos=main_window.ui.refYPos.value()
        y_width=main_window.ui.refYWidth.value()
        bg_pos=main_window.ui.bgCenter.value()
        bg_width=main_window.ui.bgWidth.value()

        if preserve_lim:
            xview=main_window.ui.x_project.canvas.ax.axis()
            yview=main_window.ui.y_project.canvas.ax.axis()
        xxlim=(0, len(xproj)-1)
        xylim=(xproj[xproj>0].min(), xproj.max()*2)
        yxlim=(0, len(yproj)-1)
        yylim=(yproj[yproj>0].min(), yproj.max()*2)

        if self._x_projection is None or len(self._x_projection.get_xdata())!=xproj.shape[0]:
            main_window.ui.x_project.clear_fig()
            main_window.ui.y_project.clear_fig()
            x_projection=main_window.ui.x_project.plot(xproj, color='blue')[0]
            main_window.ui.x_project.set_xlabel(u'x [pix]')
            main_window.ui.x_project.set_ylabel(u'I$_{max}$')
            xpos=main_window.ui.x_project.canvas.ax.axvline(x_peak, color='black')
            xleft=main_window.ui.x_project.canvas.ax.axvline(x_peak-x_width/2., color='red')
            xright=main_window.ui.x_project.canvas.ax.axvline(x_peak+x_width/2., color='red')
            xleft_bg=main_window.ui.x_project.canvas.ax.axvline(bg_pos-bg_width/2., color='black')
            xright_bg=main_window.ui.x_project.canvas.ax.axvline(bg_pos+bg_width/2., color='black')
            
            self._y_projection=main_window.ui.y_project.plot(yproj, color='blue')[0]
            main_window.ui.y_project.set_xlabel(u'y [pix]')
            main_window.ui.y_project.set_ylabel(u'I$_{max}$')
            yreg_left=main_window.ui.y_project.canvas.ax.axvline(y_pos-y_width/2., color='green')
            yreg_right=main_window.ui.y_project.canvas.ax.axvline(y_pos+y_width/2., color='green')
            ybg=main_window.ui.y_project.canvas.ax.axhline(self.y_bg, color='black')
            self.proj_lines=(xleft, xpos, xright, xleft_bg, xright_bg, yreg_left, yreg_right, ybg)
            self._x_projection = x_projection
        else:
            self._x_projection.set_ydata(xproj)
            self._y_projection.set_ydata(yproj)
            lines=self.proj_lines
            lines[0].set_xdata([x_peak-x_width/2., x_peak-x_width/2.])
            lines[1].set_xdata([x_peak, x_peak])
            lines[2].set_xdata([x_peak+x_width/2., x_peak+x_width/2.])
            lines[3].set_xdata([bg_pos-bg_width/2., bg_pos-bg_width/2.])
            lines[4].set_xdata([bg_pos+bg_width/2., bg_pos+bg_width/2.])
            lines[5].set_xdata([y_pos-y_width/2., y_pos-y_width/2.])
            lines[6].set_xdata([y_pos+y_width/2., y_pos+y_width/2.])
            lines[7].set_ydata([self.y_bg, self.y_bg])

        if preserve_lim:
            main_window.ui.x_project.canvas.ax.axis(xview)
            main_window.ui.y_project.canvas.ax.axis(yview)
        if main_window.ui.logarithmic_y.isChecked():
            main_window.ui.x_project.set_yscale('log')
            main_window.ui.y_project.set_yscale('log')
        else:
            main_window.ui.x_project.set_yscale('linear')
            main_window.ui.y_project.set_yscale('linear')
        main_window.ui.x_project.canvas.ax.set_xlim(*xxlim)
        main_window.ui.x_project.canvas.ax.set_ylim(*xylim)
        main_window.ui.y_project.canvas.ax.set_xlim(*yxlim)
        main_window.ui.y_project.canvas.ax.set_ylim(*yylim)

        main_window.ui.x_project.draw()
        main_window.ui.y_project.draw()

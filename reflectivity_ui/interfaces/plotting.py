# -*- coding: utf-8 -*-
import numpy as np

def plot_overview(self):
    '''
    X vs. Y and X vs. Tof for main channel.
    '''
    data=self._active_channel
    if data.total_counts==0:
        self.ui.xy_overview.clear()
        self.ui.xtof_overview.clear()
        
        self.ui.xy_overview.draw()
        self.ui.xtof_overview.draw()
        return
    xy=data.xydata
    xtof=data.xtofdata
    ref_norm=self.getNorm()
    if self.ui.normalizeXTof.isChecked() and ref_norm is not None:
        ref_norm=ref_norm.Rraw
        # normalize ToF dataset for wavelength distribution
        ref_norm=np.where(ref_norm>0., ref_norm, 1.)
        xtof=xtof.astype(float)/ref_norm[np.newaxis, :]
    xy_imin=xy[xy>0].min()
    xy_imax=xy.max()
    tof_imin=xtof[xtof>0].min()
    tof_imax=xtof.max()

    # for lines of the current extraction area
    x_peak=self.ui.refXPos.value()
    x_width=self.ui.refXWidth.value()
    y_pos=self.ui.refYPos.value()
    y_width=self.ui.refYWidth.value()
    bg_pos=self.ui.bgCenter.value()
    bg_width=self.ui.bgWidth.value()

    # XY plot
    if self.ui.tthPhi.isChecked():
        rad_per_pixel=data.det_size_x/data.dist_sam_det/data.xydata.shape[1]
        phi_range=xy.shape[0]*rad_per_pixel*180./np.pi
        tth_range=xy.shape[1]*rad_per_pixel*180./np.pi
        phi0=self.ui.refYPos.value()*rad_per_pixel*180./np.pi
        tth0=(data.dangle-data.dangle0)-(xy.shape[1]-data.dpix)*rad_per_pixel*180./np.pi
        self.ui.xy_overview.clear()
        if self.overview_lines is None:
            self.overview_lines=[]
        else:
            self.overview_lines=self.overview_lines[-2:]

        self.ui.xy_overview.imshow(xy, log=self.ui.logarithmic_colorscale.isChecked(),
                                 aspect='auto', cmap=self.color, origin='lower',
                                 extent=[tth_range+tth0, tth0, phi0, phi0-phi_range])
        self.ui.xy_overview.set_xlabel(u'2$\\Theta{}$ [°]')
        self.ui.xy_overview.set_ylabel(u'$\\phi{}$ [°]')
        self.ui.xy_overview.cplot.set_clim([xy_imin, xy_imax])
    else:
        self.ui.xy_overview.imshow(xy, log=self.ui.logarithmic_colorscale.isChecked(),
                                 aspect='auto', cmap=self.color, origin='lower')
        self.ui.xy_overview.set_xlabel(u'x [pix]')
        self.ui.xy_overview.set_ylabel(u'y [pix]')
        self.ui.xy_overview.cplot.set_clim([xy_imin, xy_imax])
        
        if self.overview_lines is None or len(self.overview_lines)==2:
            x1=self.ui.xy_overview.canvas.ax.axvline(x_peak-x_width/2., color='#aa0000')
            x2=self.ui.xy_overview.canvas.ax.axvline(x_peak+x_width/2., color='#aa0000')
            y1=self.ui.xy_overview.canvas.ax.axhline(y_pos-y_width/2., color='#00aa00')
            y2=self.ui.xy_overview.canvas.ax.axhline(y_pos+y_width/2., color='#00aa00')
            if self.overview_lines is not None:
                self.overview_lines=[x1, x2, y1, y2]+self.overview_lines
            else:
                self.overview_lines=[x1, x2, y1, y2]
        else:
            self.overview_lines[0].set_xdata([x_peak-x_width/2., x_peak-x_width/2.])
            self.overview_lines[1].set_xdata([x_peak+x_width/2., x_peak+x_width/2.])
            self.overview_lines[2].set_ydata([y_pos-y_width/2., y_pos-y_width/2.])
            self.overview_lines[3].set_ydata([y_pos+y_width/2., y_pos+y_width/2.])
    # XToF plot
    if self.ui.xLamda.isChecked():
        self.ui.xtof_overview.imshow(xtof[::-1], log=self.ui.logarithmic_colorscale.isChecked(),
                                     aspect='auto', cmap=self.color,
                                     extent=[data.lamda[0], data.lamda[-1], 0, data.x.shape[0]-1])
        self.ui.xtof_overview.set_xlabel(u'$\\lambda{}$ [Å]')
    else:
        self.ui.xtof_overview.imshow(xtof[::-1], log=self.ui.logarithmic_colorscale.isChecked(),
                                     aspect='auto', cmap=self.color,
                                     extent=[data.tof[0]*1e-3, data.tof[-1]*1e-3, 0, data.x.shape[0]-1])
        self.ui.xtof_overview.set_xlabel(u'ToF [ms]')
    self.ui.xtof_overview.set_ylabel(u'x [pix]')
    if len(self.overview_lines) in [0, 4]:
        x3=self.ui.xtof_overview.canvas.ax.axhline(x_peak-x_width/2., color='#aa0000')
        x4=self.ui.xtof_overview.canvas.ax.axhline(x_peak+x_width/2., color='#aa0000')
        x5=self.ui.xtof_overview.canvas.ax.axhline(bg_pos-bg_width/2., color='black')
        x6=self.ui.xtof_overview.canvas.ax.axhline(bg_pos+bg_width/2., color='black')
        self.overview_lines+=[x3, x4, x5, x6]
    else:
        self.overview_lines[-4].set_ydata([x_peak-x_width/2., x_peak-x_width/2.])
        self.overview_lines[-3].set_ydata([x_peak+x_width/2., x_peak+x_width/2.])
        self.overview_lines[-2].set_ydata([bg_pos-bg_width/2., bg_pos-bg_width/2.])
        self.overview_lines[-1].set_ydata([bg_pos+bg_width/2., bg_pos+bg_width/2.])
    self.ui.xtof_overview.cplot.set_clim([tof_imin, tof_imax])
    
    if self.ui.show_colorbars.isChecked() and self.ui.xy_overview.cbar is None:
        self.ui.xy_overview.cbar=self.ui.xy_overview.canvas.fig.colorbar(self.ui.xy_overview.cplot)
        self.ui.xtof_overview.cbar=self.ui.xtof_overview.canvas.fig.colorbar(self.ui.xtof_overview.cplot)
    self.ui.xy_overview.draw()
    self.ui.xtof_overview.draw()


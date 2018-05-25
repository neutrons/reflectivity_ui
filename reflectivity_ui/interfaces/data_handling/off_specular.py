"""
    Class to execute and hold the off-specular reflectivity calculation.
"""
import logging
import numpy as np

H_OVER_M_NEUTRON = 3.956034e-7 # h/m_n [m^2/s]


class OffSpecular(object):
    """
        Compute off-specular reflectivity
    """
    d_wavelength = 0
    Qx = None
    Qz = None
    ki_z = None
    kf_z = None
    S = None
    dS = None

    def __init__(self, cross_section_data):
        """
            :param CrossSectionData cross_section_data: processed data object
        """
        self.data_set = cross_section_data

    def __call__(self, direct_beam=None):
        """
            Extract off-specular scattering from 4D dataset (x,y,ToF,I).
            Uses a window in y to filter the 4D data
            and than sums all I values for each ToF and x channel.
            Qz,Qx,kiz,kfz is calculated using the x and ToF positions
            together with the tth-bank and direct pixel values.

            :param CrossSectionData direct_beam: if given, this data will be used to normalize the output
        """
        #TODO: correct for detector sensitivity

        x_pos = self.data_set.configuration.peak_position
        x_width = self.data_set.configuration.peak_width
        y_pos = self.data_set.configuration.low_res_position
        y_width = self.data_set.configuration.low_res_width
        scale = 1./self.data_set.proton_charge * self.data_set.configuration.scaling_factor

        # Get regions in pixels as integers
        reg = [int(round(item)) for item in [x_pos-x_width/2., x_pos+x_width/2.+1,
                                             y_pos-y_width/2., y_pos+y_width/2.+1]]

        rad_per_pixel = self.data_set.det_size_x / self.data_set.dist_sam_det / self.data_set.xydata.shape[1]

        xtth = self.data_set.direct_pixel - np.arange(self.data_set.data.shape[0])[self.data_set.active_area_x[0]:
                                                                                   self.data_set.active_area_x[1]]
        pix_offset_spec = self.data_set.direct_pixel - x_pos
        delta_dangle = self.data_set.dangle - self.data_set.dangle0
        tth_spec = delta_dangle * np.pi/180. + pix_offset_spec * rad_per_pixel
        af = delta_dangle * np.pi/180. + xtth * rad_per_pixel - tth_spec/2.
        ai = np.ones_like(af) * tth_spec / 2.

        # Background
        bck = self.data_set.get_background_vs_TOF() * scale

        v_edges = self.data_set.dist_mod_det/self.data_set.tof_edges * 1e6 #m/s
        lambda_edges = H_OVER_M_NEUTRON / v_edges * 1e10 #A

        wl = (lambda_edges[:-1] + lambda_edges[1:]) / 2.
        # The resolution for lambda is digital range with equal probability
        # therefore it is the bin size divided by sqrt(12)
        self.d_wavelength = np.abs(lambda_edges[:-1] - lambda_edges[1:]) / np.sqrt(12)
        k = 2. * np.pi / wl

        # calculate reciprocal space, incident and outgoing perpendicular wave vectors
        self.Qz=k[np.newaxis, :]*(np.sin(af)+np.sin(ai))[:, np.newaxis]
        self.Qx=k[np.newaxis, :]*(np.cos(af)-np.cos(ai))[:, np.newaxis]
        self.ki_z=k[np.newaxis, :]*np.sin(ai)[:, np.newaxis]
        self.kf_z=k[np.newaxis, :]*np.sin(af)[:, np.newaxis]

        # calculate ROI intensities and normalize by number of points
        raw_multi_dim=self.data_set.data[self.data_set.active_area_x[0]:self.data_set.active_area_x[1], reg[2]:reg[3], :]
        raw = raw_multi_dim.sum(axis=1)
        d_raw = np.sqrt(raw)

        # normalize data by width in y and multiply scaling factor
        intensity = raw/(reg[3]-reg[2]) * scale
        d_intensity = d_raw/(reg[3]-reg[2]) * scale
        self.S = intensity - bck[np.newaxis, :]
        self.dS = np.sqrt(d_intensity**2+(bck**2)[np.newaxis, :])

        if direct_beam is not None:
            if not direct_beam.configuration.tof_bins == self.data_set.configuration.tof_bins:
                logging.error("Trying to normalize with a direct beam data set with different binning")

            norm_raw_multi_dim=direct_beam.data[self.data_set.active_area_x[0]:self.data_set.active_area_x[1], reg[2]:reg[3], :]
            norm_raw = norm_raw_multi_dim.sum(axis=0).sum(axis=0)
            norm_d_raw = np.sqrt(norm_raw)
            norm_raw /= (reg[3]-reg[2]) * direct_beam.proton_charge * direct_beam.configuration.scaling_factor
            norm_d_raw /= (reg[3]-reg[2]) * direct_beam.proton_charge * direct_beam.configuration.scaling_factor

            idxs=norm_raw>0.
            self.dS[:, idxs]=np.sqrt(
                         (self.dS[:, idxs]/norm_raw[idxs][np.newaxis, :])**2+
                         (self.S[:, idxs]/norm_raw[idxs][np.newaxis, :]**2*norm_d_raw[idxs][np.newaxis, :])**2
                         )
            self.S[:, idxs]/=norm_raw[idxs][np.newaxis, :]
            self.S[:, np.logical_not(idxs)]=0.
            self.dS[:, np.logical_not(idxs)]=0.

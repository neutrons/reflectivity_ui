"""
    Data processing workflow, taking results and writing them to files.
"""
#pylint: disable=bare-except, too-many-locals
from __future__ import absolute_import, division, print_function
import sys
import os
import math
import copy
import logging
import time
import numpy as np
from ..configuration import Configuration
from . import quicknxs_io, data_manipulation, off_specular, gisans


DEFAULT_OPTIONS = dict(export_specular=True,
                       export_asym=False,
                       export_gisans=False,
                       export_offspec=False,
                       export_offspec_smooth=False,
                       format_genx=False,
                       format_matlab=False,
                       format_mantid=True,
                       format_multi=True,
                       format_numpy=False,
                       format_5cols=False,
                       output_sample_size=10,
                       output_directory='',
                       output_file_template='(instrument)_{numbers}_{item}_{state}.{type}')


class ProcessingWorkflow(object):
    """
        Carry out the reduction process for a set of data runs and manages outputs
    """
    def __init__(self, data_manager, output_options=None):
        self.data_manager = data_manager
        self.output_options = output_options if output_options else DEFAULT_OPTIONS

    def execute(self, progress=None):
        """
            Process data and write output files
            :param ProgressReporter progress: reporter object
        """
        if not self.data_manager.reduction_states:
            return

        if self.output_options['export_specular']:
            if progress is not None:
                progress(10, "Computing reflectivity")
            self.specular_reflectivity()

        if self.output_options['export_offspec'] or self.output_options['export_offspec_smooth']:
            if progress is not None:
                progress(20, "Computing off-specular reflectivity")
            self.offspec(raw=self.output_options['export_offspec'],
                         binned=self.output_options['export_offspec_smooth'])

        if progress is not None:
                progress(60, "Computing GISANS")
        if self.output_options['export_gisans']:
            self.gisans(progress=progress)

        if progress is not None:
            progress(100, "Complete")

    def get_file_name(self, run_list=None, pol_state=None, data_type='dat', process_type='Specular'):
        """
            Construct a file name according to the measurement type.
            :param list run_list: list of run numbers
            :param str pol_state: name for the polarization state
            :param str data_type: file extension
            :param str process_type: descriptor for the process type
        """
        if run_list is None:
            run_list = []
        base_name = self.output_options['output_file_template'].replace('{numbers}', '+'.join(run_list))
        base_name = base_name.replace('{instrument}', self.data_manager.active_channel.configuration.instrument.instrument_name)
        base_name = base_name.replace('{item}', process_type)
        if pol_state is not None:
            base_name = base_name.replace('{state}', pol_state)
        base_name = base_name.replace('{type}', data_type)
        return os.path.join(self.output_options['output_directory'], base_name)

    def write_quicknxs(self, output_data, output_file_base, xs=None):
        """
            Write QuickNXS output reflectivity file.
            :param dict output_data: dictionary of numpy arrays
            :param str output_file_base: template for output file paths
            :param list xs: list of cross sections available in the output_data
        """
        # Get the column names
        units = output_data['units']
        cols = output_data['columns']
        col_names = [u'%s [%s]' % (cols[i], units[i]) for i in range(len(cols))]

        # List of all output states we have to deal with
        if xs is not None:
            output_states = xs
        else:
            output_states = copy.copy(self.data_manager.reduction_states)
            if self.output_options['export_asym'] and 'SA' in output_data:
                output_states.append("SA")

        # Sanity check
        if len(output_states) == 0:
            return

        # Write out the cross-section data
        five_cols = self.output_options['format_5cols']
        for pol_state in output_states:
            # The cross-sections might have different names
            if pol_state in self.data_manager.reduction_list[0].cross_sections:
                _pol_state = self.data_manager.reduction_list[0].cross_sections[pol_state].cross_section_label
            else:
                _pol_state = pol_state
            # We might not have data for a given cross-section
            if pol_state not in output_data:
                continue

            state_output_path = output_file_base.replace('{state}', pol_state)
            quicknxs_io.write_reflectivity_header(self.data_manager.reduction_list,
                                                  self.data_manager.direct_beam_list,
                                                  state_output_path, _pol_state)
            quicknxs_io.write_reflectivity_data(state_output_path, output_data[pol_state],
                                                col_names, as_5col=five_cols)

    def write_genx(self, output_data, output_path):
        '''
            Create a Genx .gx model file with the right polarization states
            and the reflectivity data already included for convenience.
        '''
        try:
            import zipfile
            import cPickle
        except:
            logging.error("Problem importing modules: %s", sys.exc_info([1]))
            return

        # Load templates
        module_dir, _ = os.path.split(__file__)
        template_dir = os.path.join(module_dir, 'genx_templates')
        logging.error("GENX: %s", template_dir)
        if len(self.data_manager.reduction_states) == 1:
            template_path = os.path.join(template_dir, 'unpolarized.gx')
        if len(self.data_manager.reduction_states) == 2:
            template_path = os.path.join(template_dir, 'polarized.gx')
        else:
            template_path = os.path.join(template_dir, 'spinflip.gx')

        zip_template = zipfile.ZipFile(template_path, 'r')
        zip_output = zipfile.ZipFile(output_path, 'w', zip_template.compression)
        for key in ['script', 'parameters', 'fomfunction', 'config', 'optimizer']:
            zip_output.writestr(key, zip_template.read(key))

        model_data = cPickle.loads(zip_template.read('data'))
        for i, channel in enumerate(self.data_manager.reduction_states):
            if channel not in output_data:
                logging.error("Cross-section %s not in %s", channel, str(output_data.keys()))
                continue
            model_data[i].x_raw = output_data[channel][:, 0]
            model_data[i].y_raw = output_data[channel][:, 1]
            model_data[i].error_raw = output_data[channel][:, 2]
            model_data[i].xerror_raw = output_data[channel][:, 3]
            model_data[i].name = output_data['cross_sections'][channel]
            model_data[i].run_command()
        zip_output.writestr('data', cPickle.dumps(model_data, 0))
        zip_template.close()
        zip_output.close()

    def specular_reflectivity(self):
        """
            Retrieve the computed reflectivity and save it to file
        """
        # The reflectivity should always be up to date, so we don't need to recalculate it.
        # The following would be used to recalculate it:
        #    self.data_manager.calculate_reflectivity(specular=True)

        #self.data_manager.merge_data_sets(asymmetry=self.output_options['export_asym'])

        run_list = [str(item.number) for item in self.data_manager.reduction_list]

        output_data = self.get_output_data()

        # QuickNXS format
        if self.output_options['format_multi']:
            output_file_base = self.get_file_name(run_list)
            self.write_quicknxs(output_data, output_file_base)

        # Numpy arrays
        if self.output_options['format_numpy']:
            output_file = self.get_file_name(run_list, data_type='npz', pol_state='all')
            np.savez(output_file, **output_data)

        # Matlab output
        if self.output_options['format_matlab']:
            try:
                from scipy.io import savemat
                output_file = self.get_file_name(run_list, data_type='mat', pol_state='all')
                savemat(output_file, output_data, oned_as='column')
            except:
                logging.error("Could not save in matlab format: %s", sys.exc_info([1]))

        if self.output_options['format_genx']:
            output_path = self.get_file_name(run_list, data_type='gx', pol_state='all')
            self.write_genx(output_data, output_path)

        if self.output_options['format_mantid']:
            output_file = self.get_file_name(run_list, data_type='py', pol_state='all')
            script = ''
            for pol_state in self.data_manager.reduction_states:
                script += data_manipulation.generate_script(self.data_manager.reduction_list, pol_state)
            with open(output_file, 'w') as file_object:
                file_object.write(script)

    def gisans(self, progress=None):
        """
            Export GISANS.
        """
        run_list = [str(item.number) for item in self.data_manager.reduction_list]

        # Refresh the reflectivity calculation
        if progress is not None:
            progress(65, "Reducing GISANS...")

        self.data_manager.reduce_gisans(progress=None)

        if progress is not None:
            progress(75, "Binning GISANS...")

        data_dict = self.get_gisans_data(progress=None)

        if progress is not None:
            progress(90, "Writing data")

        # QuickNXS format ['smooth' is an odd name but we keep it for backward compatibility]
        output_file_base = self.get_file_name(run_list, process_type='GISANS')
        self.write_quicknxs(data_dict, output_file_base, xs=data_dict['cross_sections'].keys())

        if progress is not None:
            progress(100, "GISANS complete")

    def offspec(self, raw=True, binned=False):
        """
            Export off-specular reflectivity.
            :param bool raw: if true, the raw results will be saved
            :param bool binned: if true, the raw results will be binned and saved
        """
        run_list = [str(item.number) for item in self.data_manager.reduction_list]

        # Refresh the reflectivity calculation
        self.data_manager.reduce_offspec()

        if raw or binned:
            output_data = self.get_offspec_data()
        # Export raw result
        if raw:
            # QuickNXS format
            output_file_base = self.get_file_name(run_list, process_type='OffSpec')
            self.write_quicknxs(output_data, output_file_base)

        # Export binned result
        self.data_manager.cached_offspec = None
        if binned:
            if self.data_manager.active_channel.configuration.apply_smoothing:
                # "Smooth" version
                try:
                    smooth_output = self.smooth_offspec(output_data)
                    output_file_base = self.get_file_name(run_list, process_type='OffSpecSmooth')
                    self.write_quicknxs(smooth_output, output_file_base)
                    self.data_manager.cached_offspec = smooth_output
                except:
                    logging.error("Problem writing smooth off-spec output: %s", sys.exc_value)

            # Binned version
            binned_data, slice_data_dict = self.get_rebinned_offspec_data()
            # QuickNXS format ['smooth' is an odd name but we keep it for backward compatibility]
            output_file_base = self.get_file_name(run_list, process_type='OffSpecBinned')
            self.write_quicknxs(binned_data, output_file_base)
            if slice_data_dict is not None and 'cross_sections' in slice_data_dict:
                output_file_base = self.get_file_name(run_list, process_type='OffSpecSlice')
                self.write_quicknxs(slice_data_dict, output_file_base, xs=slice_data_dict['cross_sections'].keys())
            if self.data_manager.cached_offspec is None:
                self.data_manager.cached_offspec = binned_data

    def get_rebinned_offspec_data(self):
        """
            Get a data dictionary ready for saving
        """
        data_dict = None
        slice_data_dict = None

        # Extract common information
        if len(self.data_manager.reduction_states) == 0:
            logging.error("List of cross-sections is empty")
            return {}

        for pol_state in self.data_manager.reduction_states:
            y_list = []
            if self.output_options['off_spec_slice']:
                y_list = self.output_options['off_spec_qz_list']
            r, dr, x, y, q_data, labels = self.data_manager.rebin_offspec(pol_state,
                                                                          axes=self.output_options['off_spec_x_axis'],
                                                                          y_list=y_list,
                                                                          use_weights=self.output_options['off_spec_err_weight'],
                                                                          n_bins_x=self.output_options['off_spec_nxbins'],
                                                                          n_bins_y=self.output_options['off_spec_nybins'])
            if data_dict is None:
                data_dict = dict(units=['1/A', '1/A', 'a.u.', 'a.u.'],
                                 columns=[labels[0], labels[1], 'I', 'dI'],
                                 cross_sections={})
                slice_data_dict = dict(units=['1/A', 'a.u.', 'a.u.'],
                                       columns=[labels[0], 'I', 'dI'],
                                       cross_sections={})

            # Create array of x-values
            x_tiled = np.tile(x, len(y))
            x_tiled = x_tiled.reshape([len(y), len(x)])

            # Create array of y-values
            y_tiled = np.tile(y, len(x))
            y_tiled = y_tiled.reshape([len(x), len(y)])
            y_tiled = y_tiled.T

            rdata = np.array([x_tiled, y_tiled, r, dr]).transpose((1, 2, 0))

            if pol_state in self.data_manager.reduction_list[0].cross_sections:
                _pol_state = self.data_manager.reduction_list[0].cross_sections[pol_state].cross_section_label
            else:
                _pol_state = pol_state
            data_dict[pol_state] = [np.nan_to_num(rdata)]
            data_dict["cross_sections"][pol_state] = _pol_state

            if q_data is not None and len(q_data) > 0:
                for item in q_data:
                    slice_data_dict["cross_sections"][item[1]] = item[1]
                    slice_data_dict[item[1]] = item[0]

        return data_dict, slice_data_dict

    def get_gisans_data(self, progress=None):
        wl_npts = self.data_manager.active_channel.configuration.gisans_wl_npts
        wl_min = self.data_manager.active_channel.configuration.gisans_wl_min
        wl_max = self.data_manager.active_channel.configuration.gisans_wl_max
        qy_npts = self.data_manager.active_channel.configuration.gisans_qy_npts
        qz_npts = self.data_manager.active_channel.configuration.gisans_qz_npts
        use_pf = self.data_manager.active_channel.configuration.gisans_use_pf

        data_dict = dict(units=['1/A', '1/A', 'a.u.', 'a.u.'], cross_sections={})
        if use_pf:
            data_dict['columns'] = ['Qy', 'pf', 'I', 'dI']
        else:
            data_dict['columns'] = ['Qy', 'Qz', 'I', 'dI']

        # Extract common information
        if not self.data_manager.reduction_states or not self.data_manager.reduction_list:
            logging.error("List of cross-sections is empty")
            return data_dict

        t_0 = time.time()
        _parallel = False
        for pol_state in self.data_manager.reduction_states:
            if _parallel:
                binned_data = gisans.rebin_parallel(self.data_manager.reduction_list, pol_state,
                                                    wl_min=wl_min, wl_max=wl_max, wl_npts=wl_npts,
                                                    qy_npts=qy_npts, qz_npts=qz_npts, use_pf=use_pf)
            for i in range(wl_npts):
                wl_step = (wl_max - wl_min) / wl_npts
                _wl_min = wl_min + i * wl_step
                _wl_max = wl_min + (i + 1) * wl_step
                if _parallel:
                    _intensity, _qy, _qz_axis, _intensity_err = binned_data[i]
                else:
                    _intensity, _qy, _qz_axis, _intensity_err = self.data_manager.rebin_gisans(pol_state,
                                                                                               wl_min=_wl_min,
                                                                                               wl_max=_wl_max,
                                                                                               qy_npts=qy_npts,
                                                                                               qz_npts=qz_npts,
                                                                                               use_pf=use_pf)

                qz, qy = np.meshgrid(_qz_axis, _qy)
                rdata = np.array([qy, qz, _intensity, _intensity_err]).transpose((1, 2, 0))

                if pol_state in self.data_manager.reduction_list[0].cross_sections:
                    _pol_state = self.data_manager.reduction_list[0].cross_sections[pol_state].cross_section_label
                else:
                    _pol_state = pol_state

                _pol_state = '%.3f-%.3f_%s' % (_wl_min, _wl_max, _pol_state)
                data_dict[_pol_state] = [np.nan_to_num(rdata)]
                data_dict["cross_sections"][_pol_state] = _pol_state

        logging.info("GISANS processing time: %s sec", (time.time()-t_0))
        return data_dict

    def get_offspec_data(self):
        """
            Get a data dictionary ready for saving
        """
        data_dict = dict(units=['1/A', '1/A', '1/A', '1/A', '1/A', 'a.u.', 'a.u.'],
                         columns=['Qx', 'Qz', 'ki_z', 'kf_z', 'ki_z-kf_z', 'I', 'dI'],
                         cross_sections={})

        # Extract common information
        if not self.data_manager.reduction_states or not self.data_manager.reduction_list:
            logging.error("List of cross-sections is empty")
            return data_dict

        first_state = self.data_manager.reduction_states[0]
        p_0 = [item.cross_sections[first_state].configuration.cut_first_n_points for item in self.data_manager.reduction_list]
        p_n = [item.cross_sections[first_state].configuration.cut_last_n_points for item in self.data_manager.reduction_list]

        ki_max = 0.01
        for pol_state in self.data_manager.reduction_states:
            # The scaling factors should have been determined at this point. Just use them
            # to merge the different runs in a set.

            combined_data = []

            for item in self.data_manager.reduction_list:
                offspec = item.cross_sections[pol_state].off_spec
                Qx, Qz, ki_z, kf_z, S, dS = (offspec.Qx, offspec.Qz, offspec.ki_z, offspec.kf_z,
                                             offspec.S, offspec.dS)

                n_total = len(S[0])
                # P_0 and P_N are the number of points to cut in TOF on each side
                p_0 = item.cross_sections[pol_state].configuration.cut_first_n_points
                p_n = n_total-item.cross_sections[pol_state].configuration.cut_last_n_points

                rdata = np.array([Qx[:, p_0:p_n], Qz[:, p_0:p_n], ki_z[:, p_0:p_n], kf_z[:, p_0:p_n],
                                  ki_z[:, p_0:p_n]-kf_z[:, p_0:p_n], S[:, p_0:p_n], dS[:, p_0:p_n]]).transpose((1, 2, 0))
                combined_data.append(rdata)
                ki_max = max(ki_max, ki_z.max())

            if pol_state in self.data_manager.reduction_list[0].cross_sections:
                _pol_state = self.data_manager.reduction_list[0].cross_sections[pol_state].cross_section_label
            else:
                _pol_state = pol_state
            data_dict[pol_state] = combined_data
            data_dict["cross_sections"][pol_state] = _pol_state
        data_dict['ki_max'] = ki_max
        return data_dict

    def smooth_offspec(self, data_dict):
        """
            NOTE: 

            Create a smoothed dataset from the off-specular scattering.
            :param dict data_dict: the output of get_offspec_data()

            Note for my own integrity (MD):
               I don't think one should smooth data distributions and do any quantitative
               work with it following this process. The way this was implemented in the
               previouse QuickNXS, replicated here, is equivalent to adding an extra resolution,
               which then would have to be properly taken into account when fitting.
               In addition, the process doesn't produce errors in intensity.
               It effectively only produces a pretty picture and should only be used as such.
        """
        axes = self.data_manager.active_channel.configuration.off_spec_x_axis
        output_data=dict(cross_sections=dict())
        for channel in data_dict['cross_sections'].keys():
            data = np.hstack(data_dict[channel])
            I = data[:, :, 5].flatten()
            Qzmax = data[:, :, 2].max() * 2.

            if axes == Configuration.QX_VS_QZ:
                x=data[:, :, 0].flatten()
                y=data[:, :, 1].flatten()
                output_data['units'] = ['1/A', '1/A', 'a.u.']
                output_data['columns'] = ['Qx', 'Qz', 'I']
                axis_sigma_scaling = 2
                xysigma0 = Qzmax / 3.
            elif axes == Configuration.KZI_VS_KZF:
                x=data[:, :, 2].flatten()
                y=data[:, :, 3].flatten()
                output_data['units'] = ['1/A', '1/A', 'a.u.']
                output_data['columns'] = ['ki_z', 'kf_z', 'I']
                axis_sigma_scaling = 3
                xysigma0 = Qzmax / 6.
            else:
                x=data[:, :, 4].flatten()
                y=data[:, :, 1].flatten()
                output_data['units'] = ['1/A', '1/A', 'a.u.']
                output_data['columns'] = ['ki_z-kf_z', 'Qz', 'I']
                axis_sigma_scaling = 2
                xysigma0 = Qzmax / 3.

            x, y, I = off_specular.smooth_data(x, y, I,
                                               sigmas=self.data_manager.active_channel.configuration.off_spec_sigmas,
                                               gridx=self.data_manager.active_channel.configuration.off_spec_nxbins,
                                               gridy=self.data_manager.active_channel.configuration.off_spec_nybins,
                                               sigmax=self.data_manager.active_channel.configuration.off_spec_sigmax,
                                               sigmay=self.data_manager.active_channel.configuration.off_spec_sigmay,
                                               x1=self.data_manager.active_channel.configuration.off_spec_x_min,
                                               x2=self.data_manager.active_channel.configuration.off_spec_x_max,
                                               y1=self.data_manager.active_channel.configuration.off_spec_y_min,
                                               y2=self.data_manager.active_channel.configuration.off_spec_y_max,
                                               axis_sigma_scaling=axis_sigma_scaling, xysigma0=xysigma0)
            output_data[channel] = [np.array([x, y, I]).transpose((1, 2, 0))]
            output_data['cross_sections'][channel] = data_dict['cross_sections'][channel]
        output_data['ki_max'] = data_dict['ki_max']
        return output_data

    def get_output_data(self):
        """
            The QuickNXS format cannot be written from the merged reflectivity, so we
            have to treat it differently and give it the workspaces for each angle.
        """
        data_dict = dict(units=['1/A', 'a.u.', 'a.u.', '1/A', 'rad'],
                         columns=['Qz', 'R', 'dR', 'dQz', 'theta'],
                         cross_sections={})

        # Extract common information
        if not self.data_manager.reduction_states or not self.data_manager.reduction_list:
            logging.error("List of cross-sections is empty")
            return data_dict
        first_state = self.data_manager.reduction_states[0]
        p_0 = [item.cross_sections[first_state].configuration.cut_first_n_points for item in self.data_manager.reduction_list]
        p_n = [item.cross_sections[first_state].configuration.cut_last_n_points for item in self.data_manager.reduction_list]

        for pol_state in self.data_manager.reduction_states:
            # The scaling factors should have been determined at this point. Just use them
            # to merge the different runs in a set.
            ws_list = data_manipulation.get_scaled_workspaces(self.data_manager.reduction_list, pol_state)

            # If the reflectivity calculation failed, we may not have data to work with
            # for this cross-section.
            if len(ws_list) == 0:
                continue

            combined_data = []
            for i, ws in enumerate(ws_list):
                _x = ws.readX(0)
                n_total = len(_x)
                x = ws.readX(0)[p_0[i]:n_total-p_n[i]]
                y = ws.readY(0)[p_0[i]:n_total-p_n[i]]
                dy = ws.readE(0)[p_0[i]:n_total-p_n[i]]
                dx = ws.readDx(0)[p_0[i]:n_total-p_n[i]]
                tth_value = ws.getRun().getProperty("SANGLE").getStatistics().mean * math.pi / 180.0
                tth = np.ones(len(x)) * tth_value
                combined_data.append(np.vstack((x, y, dy, dx, tth)).transpose())

            _output_data = np.vstack(combined_data)
            ordered = np.argsort(_output_data, axis=0).transpose()[0]
            output_data = _output_data[ordered]

            if pol_state in self.data_manager.reduction_list[0].cross_sections:
                _pol_state = self.data_manager.reduction_list[0].cross_sections[pol_state].cross_section_label
            else:
                _pol_state = pol_state
            data_dict[pol_state] = output_data
            data_dict["cross_sections"][pol_state] = _pol_state

        # Asymmetry
        if self.output_options['export_asym']:
            p_state, m_state = self.data_manager.determine_asymmetry_states()
            if p_state and m_state:
                # Get the list of workspaces
                asym_data = []
                if p_state in data_dict and m_state in data_dict:
                    # Sometimes the number of points may be different if the last few points had no signal.
                    for i in range(len(data_dict[p_state])):
                        p_point = data_dict[p_state][i]
                        i_m = len(data_dict[m_state][data_dict[m_state].T[0] < p_point[0]])
                        m_point = data_dict[m_state][i_m]
                        if p_point[1] > 0 and m_point[1] > 0:
                            ratio = (p_point[1] - m_point[1]) / (p_point[1] + m_point[1])
                            d_ratio = 2.0 / (p_point[1] + m_point[1])**2
                            d_ratio *= math.sqrt(m_point[1]**2 * p_point[2]**2 + p_point[1]**2 * m_point[2]**2)
                            asym_data.append([p_point[0], ratio, d_ratio, p_point[3], p_point[4]])
                    data_dict['SA'] = np.asarray(asym_data)
                else:
                    logging.error("Asym request but failed: %s %s %s %s", p_state, m_state,
                                  len(data_dict[p_state]), len(data_dict[m_state]))

        return data_dict

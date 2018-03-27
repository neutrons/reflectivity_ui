"""
    Data processing workflow, taking results and writing them to files.
    #TODO: write mantid script
"""
#pylint: disable=bare-except
from __future__ import absolute_import, division, print_function
import sys
import os
import math
import copy
import numpy as np
import logging
from . import quicknxs_io, data_manipulation

# Standard names for array outputs
STD_CHANNELS={'x': 'unpolarized',
              '+': 'up',
              '-': 'down',
              '++': 'upup',
              '--': 'downdown',
              '+-': 'updown',
              '-+': 'downup',
              }


class ProcessingWorkflow(object):
    """
        Carry out the reduction process for a set of data runs and manages outputs
    """
    def __init__(self, data_manager, output_options):
        self.data_manager = data_manager
        self.output_options = output_options
        self.output_data_dict = {}

    def execute(self, progress=None):
        """
            Process data and write output files
            :param ProgressReporter progress: reporter object
        """
        if len(self.data_manager.reduction_states) == 0:
            return

        progress.create()
        if self.output_options['export_specular']:
            if progress is not None:
                progress(10, "Computing reflectivity")
            self.specular_reflectivity()

        if self.output_options['export_offspec']:
            pass

        if self.output_options['export_offspec_corr']:
            pass

        if self.output_options['export_offspec_smooth']:
            pass

        if self.output_options['export_gisans']:
            pass
        progress(100, "Complete")

    def get_file_name(self, run_list=[], pol_state=None, data_type='dat', process_type='Specular'):
        """
            Construct a file name according to the measurement type.
            :param list run_list: list of run numbers
            :param str pol_state: name for the polarization state
            :param str data_type: file extension
            :param str process_type: descriptor for the process type
        """
        base_name = self.output_options['output_file_template'].replace('{numbers}', '+'.join(run_list))
        base_name = base_name.replace('{instrument}', self.data_manager.active_channel.configuration.instrument.instrument_name)
        base_name = base_name.replace('{item}', process_type)
        if pol_state is not None:
            base_name = base_name.replace('{state}', pol_state)
        base_name = base_name.replace('{type}', data_type)
        return os.path.join(self.output_options['output_directory'], base_name)

    def write_quicknxs(self, output_data, output_file_base):
        """
            Write QuickNXS output reflectivity file.
            :param dict output_data: dictionary of numpy arrays
            :param str output_file_base: template for output file paths
        """
        # List of all output states we have to deal with
        output_states = copy.copy(self.data_manager.reduction_states)
        if self.output_options['export_asym'] and 'SA' in output_data:
            output_states.append("SA")

        # Create combined file header
        combined_output_path = output_file_base.replace('{state}', 'all')
        if self.output_options['format_combined']:
            all_states = ', '.join(self.data_manager.reduction_states)
            quicknxs_io.write_reflectivity_header(self.data_manager.reduction_list, combined_output_path, all_states)

        # Write out the cross-section data
        for pol_state in output_states:

            # The cross-sections might have different names
            output_xs_name = STD_CHANNELS.get(pol_state, pol_state)
            # We might not have data for a given cross-section
            if output_xs_name not in output_data:
                continue

            if self.output_options['format_combined']:
                quicknxs_io.write_reflectivity_data(combined_output_path, output_data[output_xs_name], pol_state, as_multi=True)

            if self.output_options['format_multi']:
                state_output_path = output_file_base.replace('{state}', pol_state) #self.get_file_name(run_list, pol_state=pol_state)
                quicknxs_io.write_reflectivity_header(self.data_manager.reduction_list, state_output_path, pol_state)
                quicknxs_io.write_reflectivity_data(state_output_path, output_data[output_xs_name], pol_state, as_multi=False)

    def write_genx(self, output_data, output_file_base):
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

        output_path = output_file_base.replace('{state}', 'all')

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

        iz=zipfile.ZipFile(template_path, 'r')
        oz=zipfile.ZipFile(output_path, 'w', iz.compression)
        for key in ['script', 'parameters', 'fomfunction', 'config', 'optimizer']:
            oz.writestr(key, iz.read(key))

        model_data=cPickle.loads(iz.read('data'))
        for i, channel in enumerate(self.data_manager.reduction_states):
            model_data[i].x_raw=output_data[channel][:, 0]
            model_data[i].y_raw=output_data[channel][:, 1]
            model_data[i].error_raw=output_data[channel][:, 2]
            model_data[i].xerror_raw=output_data[channel][:, 3]
            model_data[i].name=channel
            model_data[i].run_command()
        oz.writestr('data', cPickle.dumps(model_data, 0)) # dup as version 2 pickle
        iz.close()
        oz.close()

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
        output_file_base =  self.get_file_name(run_list)
        self.write_quicknxs(output_data, output_file_base)

        # Numpy arrays
        if self.output_options['format_numpy']:
            output_file =  self.get_file_name(run_list, data_type='npz', pol_state='all')
            np.savez(output_file, **output_data)

        if self.output_options['format_plot']: pass

        # Matlab output
        if self.output_options['format_matlab']:
            try:
                from scipy.io import savemat
                output_file =  self.get_file_name(run_list, data_type='mat', pol_state='all')
                savemat(output_file, output_data, oned_as='column')
            except:
                logging.error("Could not save in matlab format: %s", sys.exc_info([1]))

        if self.output_options['format_genx']:
            self.write_genx(output_data, output_file_base)

    def get_output_data(self):
        """
            The QuickNXS format cannot be written from the merged reflectivity, so we
            have to treat it differently and give it the workspaces for each angle.
        """
        data_dict = dict(units=['1/A', 'a.u.', 'a.u.', '1/A', 'rad'],
                         columns=['Qz', 'R', 'dR', 'dQz', 'theta'])

        # Extract common information
        if len(self.data_manager.reduction_states) == 0:
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
            for i in range(len(ws_list)):
                ws = ws_list[i]
                _x = ws.readX(0)
                n_total = len(_x)
                x = ws.readX(0)[p_0[i]:n_total-p_n[i]]
                y = ws.readY(0)[p_0[i]:n_total-p_n[i]]
                dy = ws.readE(0)[p_0[i]:n_total-p_n[i]]
                dx = ws.readDx(0)[p_0[i]:n_total-p_n[i]]
                tth_value = ws.getRun().getProperty("SANGLE").getStatistics().mean * math.pi / 180.0
                tth = np.ones(len(x)) * tth_value
                combined_data.append(np.vstack((x,y,dy,dx,tth)).transpose())

            _output_data = np.vstack(combined_data)
            ordered = np.argsort(_output_data, axis=0).transpose()[0]
            output_data = _output_data[ordered]

            output_xs_name = STD_CHANNELS.get(pol_state, pol_state)
            data_dict[output_xs_name] = output_data

        # Asymmetry
        if self.output_options['export_asym']:
            p_state, m_state = self.data_manager.determine_asymmetry_states()
            p_state = STD_CHANNELS.get(p_state, p_state)
            m_state = STD_CHANNELS.get(m_state, m_state)

            # Get the list of workspaces
            asym_data = []
            if p_state in data_dict and m_state in data_dict \
                and len(data_dict[p_state]) == len(data_dict[m_state]):

                for i in range(len(data_dict[p_state])):
                    p_point = data_dict[p_state][i]
                    m_point = data_dict[m_state][i]

                    ratio = (p_point[1] - m_point[1]) / (p_point[1] + m_point[1])
                    d_ratio = 2.0 / (p_point[1] + m_point[1])**2
                    d_ratio *= math.sqrt(m_point[1]**2 * p_point[2]**2 + p_point[1]**2 * m_point[2]**2)

                    asym_data.append([p_point[0], ratio, d_ratio, p_point[3], p_point[4]])
    
                data_dict['SA'] = np.asarray(asym_data)
            else:
                logging.error("Asym request but failed: %s %s %s %s", p_state, m_state,
                              len(data_dict[p_state]), len(data_dict[m_state]))

        return data_dict
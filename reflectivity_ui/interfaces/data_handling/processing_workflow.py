"""
    Data processing workflow, taking results and writing them to files.
    #TODO: write mantid script
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
from . import quicknxs_io, data_manipulation

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

    def asymmetry(self):
        """
            Determine which cross-section to use to compute asymmetry, compute it, and write it to file.
        """
        # Inspect cross-section
        # - For two states, just calculate the asymmetry using those two
        p_state = None
        m_state = None
        if len(self.data_manager.reduction_states) == 2:
            p_state = self.data_manager.reduction_states[0]
            m_state = self.data_manager.reduction_states[1]

        # - For the traditional four states, pick the right ones by hand
        elif len(self.data_manager.reduction_states) == 4:
            if '++' in self.data_manager.reduction_states \
            and '--' in self.data_manager.reduction_states:
                p_state = '++'
                m_state = '--'

        # - If we haven't made sense of it yet, take the first and last cross-sections 
        if p_state is None and m_state is None and len(self.data_manager.reduction_states)>=2:
            p_state = self.data_manager.reduction_states[0]
            m_state = self.data_manager.reduction_states[-1]

        if p_state is None or m_state is None:
            return

        # Get the list of workspaces
        ws_list = []
        for item in self.data_manager.reduction_list:
            p_ws = item.cross_sections[p_state].reflectivity_workspace
            m_ws = item.cross_sections[m_state].reflectivity_workspace
            ratio_ws = (p_ws - m_ws) / (p_ws + m_ws)
            ws_list.append(ratio_ws)

        self._write_quicknxs(ws_list, pol_state='SA', create_combined=False)

    def get_file_name(self, run_list=[], pol_state='', data_type='dat', process_type='Specular'):
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
        base_name = base_name.replace('{state}', pol_state)
        base_name = base_name.replace('{type}', data_type)
        return os.path.join(self.output_options['output_directory'], base_name)

    def write_quicknxs(self):
        """
            The QuickNXS format cannot be written from the merged reflectivity, so we
            have to treat it differently and give it the workspaces for each angle.
        """
        combined_file_created = False
        for pol_state in self.data_manager.reduction_states:
            # The scaling factors should have been determined at this point. Just use them
            # to merge the different runs in a set.
            ws_list = data_manipulation.get_scaled_workspaces(self.data_manager.reduction_list, pol_state)

            trim_first = [item.cross_sections[pol_state].configuration.cut_first_n_points for item in self.data_manager.reduction_list]
            trim_last = [item.cross_sections[pol_state].configuration.cut_last_n_points for item in self.data_manager.reduction_list]

            # Append data to QuickNXS file
            self._write_quicknxs(ws_list, pol_state=pol_state, create_combined=not combined_file_created)
            combined_file_created = True

        # Call WriteQuickNXS, with an option to compute the asymmetry from two specified x-sections
        # Specify which ones with an input paramater: ASYM_INDICES = [0,3]

        if self.output_options['export_asym']:
            self.asymmetry()

    def _write_quicknxs(self, ws_list, pol_state='', create_combined=True):
        """
            Write reflectivity in QuickNXS format
            :param list ws_list: list of reduced mantid workspaces
            :param str pol_state: name forthe polarization state
            :param bool create_combined: if True, a file for the combined states will be created
        """
        # Get the list of runs
        run_list = [str(ws.getRunNumber()) for ws in ws_list]

        # Create a name for the combined file
        combined_output_path = self.get_file_name(run_list, pol_state='all')

        # Create a name for each cross-section file
        output_path = self.get_file_name(run_list, pol_state=pol_state)

        exporter = quicknxs_io.Exporter(ws_list)
        if self.output_options['format_combined']:
            if create_combined:
                all_states = ', '.join(self.data_manager.reduction_states)
                exporter.write_reflectivity_header(combined_output_path, all_states)
            exporter.write_reflectivity_data(combined_output_path, pol_state, as_multi=True)

        if self.output_options['format_multi']:
            exporter.write_reflectivity_header(output_path, pol_state)
            exporter.write_reflectivity_data(output_path, pol_state, as_multi=False)

    def specular_reflectivity(self):
        """
            Retrieve the computed reflectivity and save it to file
        """
        # The reflectivity should always be up to date, so we don't need to recalculate it.
        # The following would be used to recalculate it:
        #    self.data_manager.calculate_reflectivity(specular=True)

        self.data_manager.merge_data_sets(asymmetry=self.output_options['export_asym'])

        # The QuickNXS format is treated differently
        self.write_quicknxs()


        if self.output_options['format_numpy']: pass

        if self.output_options['format_plot']: pass

        if self.output_options['format_matlab']: pass

        if self.output_options['format_genx']: pass


"""
    Data manager. Holds information about the current data location
    and manages the data cache.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import logging
from reflectivity_ui.interfaces.data_handling.data_set import NexusData

class DataManager(object):
    MAX_CACHE = 50

    def __init__(self, current_directory):
        self.current_directory = current_directory
        self.current_file_name = None
        # Current data set
        self._nexus_data = None
        self.active_channel = None
        # Cache of loaded data
        self._cache = []

        # The following is information about the data to be combined together
        # List of data sets
        self.reduction_list = []
        self.direct_beam_list = []
        # List of cross-sections common to all reduced data sets
        self.reduction_states = []

    @property
    def data_sets(self):
        if self._nexus_data is None:
            return None
        return self._nexus_data.cross_sections

    @property
    def current_file(self):
        if self._nexus_data is None:
            return None
        return self._nexus_data.file_path

    def get_cachesize(self):
        return sum([d.nbytes for d in self._cache])

    def clear_cache(self):
        self._cache = []

    def set_active_data_from_reduction_list(self, index):
        """
            Set a data set in the reduction list as the active
            data set according to its index.
            :param int index: index in the reduction list
        """
        if index < len(self.reduction_list):
            self._nexus_data = self.reduction_list[index]
            self.set_channel(0)

    def set_active_data_from_direct_beam_list(self, index):
        """
            Set a data set in the direct beam list as the active
            data set according to its index.
            :param int index: index in the direct beam list
        """
        if index < len(self.direct_beam_list):
            self._nexus_data = self.direct_beam_list[index]
            self.set_channel(0)

    def set_channel(self, index):
        """
            Set the current channel to the specified index, or zero
            if it doesn't exist.
        """
        channels = self.data_sets.keys()
        if index < len(channels):
            self.active_channel = self.data_sets[channels[index]]
            return True
        else:
            self.active_channel = self.data_sets[channels[0]]
        return False

    def is_active(self, data_set):
        """
            Returns True of the given data set is the active data set.
            :param NexusData: data set object
        """
        return data_set == self._nexus_data

    def is_active_data_compatible(self):
        """
            Determine whether the currently active data set is compatible
            with the data sets that are currently part of the reduction list.
        """
        return self.reduction_states == self.data_sets.keys() or self.reduction_list == []

    def find_data_in_reduction_list(self, nexus_data):
        """
            Look for the given data in the reduction list.
            Return the index within the reduction list or none.
            :param NexusData: data set object
        """
        for i in range(len(self.reduction_list)):
            if nexus_data == self.reduction_list[i]:
                return i
        return None

    def find_data_in_direct_beam_list(self, nexus_data):
        """
            Look for the given data in the direct beam list.
            Return the index within the direct beam list or none.
            :param NexusData: data set object
        """
        for i in range(len(self.direct_beam_list)):
            if nexus_data == self.direct_beam_list[i]:
                return i
        return None

    def find_active_data_id(self):
        """
            Look for the active data in the reduction list.
            Return the index within the reduction list or none.
        """
        return self.find_data_in_reduction_list(self._nexus_data)

    def find_active_direct_beam_id(self):
        """
            Look for the active data in the direct beam list.
            Return the index within the direct beam list or none.
        """
        return self.find_data_in_direct_beam_list(self._nexus_data)

    def add_active_to_reduction(self):
        """
            Add active data set to reduction list
        """
        if not self._nexus_data in self.reduction_list:
            if self.is_active_data_compatible():
                if len(self.reduction_list) == 0:
                    self.reduction_states = self.data_sets.keys()
                self.reduction_list.append(self._nexus_data)
                return True
            else:
                logging.error("The data you are trying to add has different cross-sections")
        return False

    def add_active_to_normalization(self):
        """
            Add active data set to the direct beam list
        """
        if not self._nexus_data in self.direct_beam_list:
            self.direct_beam_list.append(self._nexus_data)
            return True
        return False

    def remove_active_from_normalization(self):
        """
            Remove the active data set from the direct beam list
        """
        for i in range(len(self.direct_beam_list)):
            if self.direct_beam_list[i] == self._nexus_data:
                self.direct_beam_list.pop(i)
                return i
        return -1

    def clear_direct_beam_list(self):
        """
            Remove all items from the direct beam list, and make
            sure to remove links to those items in the scattering data sets.

            TODO: remove links from scattering data sets.
        """
        self.direct_beam_list = []

    def load(self, file_path, configuration, force=False):
        """
            Load a data file
            :param str file_path: file path
            :param Configuration configuration: configuration to use to load the data
            :param bool force: it True, existing data will be replaced if it exists.
        """
        nexus_data = None
        is_from_cache = False
        reduction_list_id = None
        direct_beam_list_id = None

        # Check whether the file is in cache
        for i in range(len(self._cache)):
            if self._cache[i].file_path == file_path:
                if force:
                    # Check whether the data is in the reduction list before
                    # removing it.
                    reduction_list_id = self.find_data_in_reduction_list(self._cache[i])
                    direct_beam_list_id = self.find_data_in_direct_beam_list(self._cache[i])
                    self._cache.pop(i)
                else:
                    nexus_data = self._cache[i]
                    is_from_cache = True
                break

        # If we don't have the data, load it
        if nexus_data is None:
            nexus_data = NexusData(file_path, configuration)
            nexus_data.load()

        if nexus_data is not None:
            self._nexus_data = nexus_data
            directory, file_name = os.path.split(file_path)
            self.current_directory = directory
            self.current_file_name = file_name

            # Find suitable direct beam
            if configuration.match_direct_beam:
                self.find_best_direct_beam()

            # If we didn't get this data set from our cache,
            # then add it and compute its reflectivity.
            if not is_from_cache:
                # Replace reduction and normalization entries as needed
                if reduction_list_id is not None:
                    self.reduction_list[reduction_list_id] = nexus_data
                if direct_beam_list_id is not None:
                    self.direct_beam_list[direct_beam_list_id] = nexus_data
                # Compute reflectivity
                nexus_data.calculate_reflectivity()
                while len(self._cache)>=self.MAX_CACHE:
                    self._cache.pop(0)
                self._cache.append(nexus_data)
            return self.data_sets
        logging.error("Nothing to load for file %s", file_path)
        return None

    def update_configuration(self, configuration, active_only=False, nexus_data=None):
        """
            Update configuration
        """
        if active_only:
            self.active_channel.update_configuration(configuration)
        elif nexus_data is not None:
            nexus_data.update_configuration(configuration)
        else:
            self._nexus_data.update_configuration(configuration)

    def calculate_reflectivity(self, configuration=None, active_only=False, nexus_data=None):
        """
            Calculater reflectivity using the current configuration
        """
        # Try to find the direct beam in the list of direct beam data sets
        direct_beam = None
        
        # Select the data to work on
        if nexus_data is None:
            nexus_data = self._nexus_data

        # Get the direct beam info from the configuration
        # All the cross sections should have the same direct beam file.
        data_keys = nexus_data.cross_sections.keys()
        if len(data_keys) == 0:
            logging.error("DataManager.calculate_reflectivity: nothing to compute")
            return

        data_xs = nexus_data.cross_sections[data_keys[0]]

        if data_xs.configuration is not None and data_xs.configuration.normalization is not None:
            for item in self.direct_beam_list:
                if item.number == data_xs.configuration.normalization:
                    keys = item.cross_sections.keys()
                    if len(keys) >= 1:
                        if len(keys) > 1:
                            logging.error("More than one cross-section for the direct beam, using the first one")
                        direct_beam = item.cross_sections[keys[0]]
            if direct_beam is None:
                logging.error("The specified direct beam is not available: skipping")

        if active_only:
            self.active_channel.reflectivity(direct_beam=direct_beam, configuration=configuration)
        else:
            nexus_data.calculate_reflectivity(direct_beam=direct_beam, configuration=configuration)

    def find_best_direct_beam(self):
        """
            Find the best direct beam in the direct beam list for the active data
            Returns a run number.
            Returns True if we have updated the data with a new normalization run.
        """
        closest = None
        for item in self.direct_beam_list:
            channel = item.cross_sections[item.cross_sections.keys()[0]]
            if self.active_channel.configuration.instrument.direct_beam_match(self.active_channel, channel):
                if closest is None:
                    closest = item.number
                elif abs(item.number-self.active_channel.number) < abs(closest-self.active_channel.number):
                    closest = item.number

        if closest is None:
            # If we didn't find a direct beam, try with just the wavelength
            for item in self.direct_beam_list:
                channel = item.cross_sections[item.cross_sections.keys()[0]]
                if self.active_channel.configuration.instrument.direct_beam_match(self.active_channel, channel, skip_slits=True):
                    if closest is None:
                        closest = item.number
                    elif abs(item.number-self.active_channel.number) < abs(closest-self.active_channel.number):
                        closest = item.number
        return self._nexus_data.set_parameter("normalization", closest)

"""
    Data manager. Holds information about the current data location
    and manages the data cache.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
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

    def is_active_data_compatible(self):
        """
            Determine whether the currently active data set is compatible
            with the data sets that are currently part of the reduction list.
        """
        return self.reduction_states == self.data_sets.keys() or self.reduction_list == []

    def add_active_to_reduction(self):
        """
            Add active data set to reduction list
        """
        if not self._nexus_data in self.reduction_list:
            if self.is_active_data_compatible():
                self.reduction_list.append(self._nexus_data)
            else:
                logging.error("The data you are trying to add has different cross-sections")
            return True
        if len(self.reduction_list) == 1:
            self.reduction_states = self.data_sets.keys()
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
        """
        nexus_data = None
        # Check whether the file is in cache
        for i in range(len(self._cache)):
            if self._cache[i].file_path == file_path:
                if force:
                    self._cache.pop(i)
                else:
                    nexus_data = self._cache[i]
                break

        # If we don't have the data, load it
        if nexus_data is None:
            nexus_data = NexusData(file_path, configuration)
            nexus_data.load()

        if nexus_data is not None:
            self._nexus_data = nexus_data
            while len(self._cache)>=self.MAX_CACHE:
                self._cache.pop(0)
            self._cache.append(nexus_data)

            directory, file_name = os.path.split(file_path)
            self.current_directory = directory
            #self.current_file = file_path
            self.current_file_name = file_name
            #self.data_sets = nexus_data.cross_sections
            return self.data_sets
        logging.error("Nothing to load for file %s", file_path)
        return None

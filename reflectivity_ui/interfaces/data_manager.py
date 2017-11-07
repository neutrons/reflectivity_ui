"""
    Data manager. Holds information about the current data location
    and manages the data cache.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import os
import logging
from .data_handling.loader import NexusData

class DataManager(object):
    MAX_CACHE = 50

    def __init__(self, current_directory):
        self.current_directory = current_directory
        self.current_file_name = None
        # Current data set
        self._nexus_data = None
        # Cache of loaded data
        self._cache = []

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

    def load(self, file_path, configuration, force=False):
        """
            Load a data file
            :param str file_path: file path
        """
        nexus_data = None
        # Only load from cache is force is false.
        if not force:
            # Check whether the file is in cache
            for item in self._cache:
                if item.file_path == file_path:
                    nexus_data = item
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

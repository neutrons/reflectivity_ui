"""
    Data manager. Holds information about the current data location
    and manages the data cache.
"""

class DataManager(object):
    def __init__(self, current_directory):
        self.current_directory = current_directory
        self.current_file = None
        self.current_file_name = None
        self.active_channel = None
        self.data_sets = None
        self._cache = []


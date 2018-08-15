import unittest
import sys
sys.path.append('..')
import os

from reflectivity_ui.interfaces.data_handling.quicknxs_io import read_reduced_file
from reflectivity_ui.interfaces.data_manager import DataManager
from reflectivity_ui.interfaces.configuration import Configuration

class DataLoaderTest(unittest.TestCase):
    def test_simple_load(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_28613+28614+28615+28616+28617+28618+28619_Specular_++.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 7)
        self.assertEqual(len(data_list), 7)
        self.assertEqual(data_list[0][2].peak_position, 179.5)
        self.assertEqual(data_list[0][2].normalization, 28610)

    def test_load_no_db(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_29160_Specular_++.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 0)
        self.assertEqual(len(data_list), 1)

    def test_load_from_ar(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_29526_Off_Off_combined_autoreduced.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 4)
        self.assertEqual(len(data_list), 4)

    def test_load_from_quicknxs(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_29526_quicknxs.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 4)
        self.assertEqual(len(data_list), 4)

    def test_load_from_mismatch(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_29782_empty_db.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 5)
        self.assertEqual(len(data_list), 6)
        self.assertEqual(data_list[4][2].normalization, None)

class DataManagerTest(unittest.TestCase):

    def test_manager(self):
        manager = DataManager(os.getcwd())
        manager.load("REF_M_29160", Configuration())

        self.assertEqual(manager.current_file, "REF_M_29160")

        manager.add_active_to_reduction()
        self.assertEqual(manager.find_data_in_reduction_list(manager._nexus_data), 0)
        self.assertEqual(manager.find_data_in_direct_beam_list(manager._nexus_data), None)

        self.assertTrue(manager.add_active_to_normalization())
        self.assertEqual(manager.remove_active_from_normalization(), 0)

        manager.set_active_data_from_reduction_list(0)
        manager.set_active_data_from_direct_beam_list(0)
        manager.calculate_reflectivity()

    def test_load_reduced(self):
        manager = DataManager(os.getcwd())
        manager.load_data_from_reduced_file('data/REF_M_29160_Specular_++.dat')

if __name__ == '__main__':
    unittest.main()

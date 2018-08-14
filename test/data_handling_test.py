import unittest
import sys
sys.path.append('..')
import os

from reflectivity_ui.interfaces.data_handling.quicknxs_io import read_reduced_file


class DataLoaderTest(unittest.TestCase):
    def test_simple_load(self):
        file_path = os.path.join(os.path.dirname(__file__), 'data',
                                 'REF_M_28613+28614+28615+28616+28617+28618+28619_Specular_++.dat')
        db_list, data_list = read_reduced_file(file_path)
        self.assertEqual(len(db_list), 7)
        self.assertEqual(len(data_list), 7)
        self.assertEqual(data_list[0][2].peak_position, 179.5)
        self.assertEqual(data_list[0][2].normalization, 28610)

if __name__ == '__main__':
    unittest.main()

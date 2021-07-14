from __future__ import absolute_import, division, print_function, unicode_literals

import os
import sys
import unittest

import reflectivity_ui.interfaces.data_handling.data_manipulation as dm
from reflectivity_ui.interfaces.data_handling.filepath import RunNumbers, FilePath
from reflectivity_ui.interfaces.data_handling import ApplicationConfiguration
from reflectivity_ui.interfaces.data_handling.quicknxs_io import read_reduced_file
from reflectivity_ui.interfaces.data_manager import DataManager
from reflectivity_ui.interfaces.configuration import Configuration

sys.path.append('..')


class ApplicationConfigurationTest(unittest.TestCase):
    def test_init(self):
        # System's installation of mantid
        application_conf = ApplicationConfiguration()
        assert os.path.dirname(application_conf.mantid_path) in sys.path
        # Custom "installation" of mantid
        mantid_path = '/tmp/mantid41'
        if not os.path.isdir(mantid_path):
            os.makedirs('/tmp/mantid41')
        application_conf = ApplicationConfiguration(root_dir='/tmp')
        assert application_conf.mantid_path == '/tmp/mantid41'
        assert application_conf.mantid_version == '4.1.0'


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

        q_range = manager._nexus_data.get_q_range()
        self.assertAlmostEqual(q_range[0], 0.034, delta=0.05)
        self.assertAlmostEqual(q_range[1], 0.068, delta=0.05)

        self.assertTrue(manager.add_active_to_normalization())
        self.assertEqual(manager.remove_active_from_normalization(), 0)

        manager.set_active_data_from_reduction_list(0)
        manager.set_active_data_from_direct_beam_list(0)
        manager.calculate_reflectivity()
        manager.calculate_reflectivity(specular=False)
        manager.strip_overlap()

        dm.generate_script(manager.reduction_list, manager.reduction_states[0])
        dm.stitch_reflectivity(manager.reduction_list)
        dm.merge_reflectivity(manager.reduction_list, manager.reduction_states[0])
        dm.get_scaled_workspaces(manager.reduction_list, manager.reduction_states[0])
        dm.stitch_reflectivity(manager.reduction_list)

    def test_load_reduced(self):
        manager = DataManager(os.getcwd())
        manager.load_data_from_reduced_file('data/REF_M_29160_Specular_++.dat')


class RunNumberTest(unittest.TestCase):
    def test_init(self):
        self.assertEquals(RunNumbers(123).numbers, [123])
        self.assertEquals(RunNumbers('123').numbers, [123])
        self.assertEquals(RunNumbers([123, '126', 125]).numbers, [123, 125, 126])
        self.assertEquals(RunNumbers('7:10+3:5+1').numbers, [1, 3, 4, 5, 7, 8, 9, 10])
        self.assertEquals(RunNumbers('7:10 + 3:5 + 1').numbers, [1, 3, 4, 5, 7, 8, 9, 10])

    def test_long(self):
        runs = RunNumbers([7, 8, 9, 10, 3, 4, 5, 1])
        self.assertEqual(runs.long, '1+3+4+5+7+8+9+10')

    def test_short(self):
        runs = RunNumbers([7, 8, 9, 10, 3, 4, 5, 1])
        self.assertEqual(runs.short, '1+3:5+7:10')


class FilePathTest(unittest.TestCase):
    def test_init(self):
        self.assertEqual(FilePath(u'/SNS/REF_M_1.nxs').path, u'/SNS/REF_M_1.nxs')
        self.assertEqual(FilePath([u'/SNS/REF_M_2.nxs', u'/SNS/REF_M_1.nxs']).path,
                         u'/SNS/REF_M_1.nxs+/SNS/REF_M_2.nxs')
        self.assertEqual(FilePath([u'/SNS/REF_M_2.nxs', u'/SNS/REF_M_1.nxs'], sort=False).path,
                         u'/SNS/REF_M_2.nxs+/SNS/REF_M_1.nxs')
        self.assertEqual(FilePath(u'/SNS/REF_M_2.nxs+/SNS/REF_M_1.nxs').path,
                         u'/SNS/REF_M_1.nxs+/SNS/REF_M_2.nxs')

    def test_join(self):
        self.assertEqual(FilePath.join(u'/SNS', u'REF_M_1.nxs'), u'/SNS/REF_M_1.nxs')
        self.assertEqual(FilePath.join(u'/SNS', u'REF_M_2.nxs+REF_M_1.nxs'),
                         u'/SNS/REF_M_1.nxs+/SNS/REF_M_2.nxs')

    def test_unique_dirname(self):
        self.assertTrue(FilePath.unique_dirname(u'/SNS/REF_M_1.nxs+/SNS/REF_M_2.nxs'))
        self.assertEqual(FilePath.unique_dirname(u'/NSN/REF_M_1.nxs+/SNS/REF_M_2.nxs'), False)

    def test_single_paths(self):
        self.assertEquals(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').single_paths,
                          [u'/SNS/REF_M_1.nxs', u'/SNS/REF_M_3.nxs'])

    def test_is_composite(self):
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs').is_composite, False)
        self.assertTrue(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').is_composite)

    def test_dirname(self):
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs').dirname, u'/SNS')
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').dirname, u'/SNS')

    def test_basename(self):
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs').basename, u'REF_M_3.nxs')
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').basename, u'REF_M_1.nxs+REF_M_3.nxs')

    def test_first_path(self):
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs').first_path, u'/SNS/REF_M_3.nxs')
        self.assertEqual(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').first_path, u'/SNS/REF_M_1.nxs')

    def test_split(self):
        self.assertEquals(FilePath(u'/SNS/REF_M_3.nxs').split(), (u'/SNS', u'REF_M_3.nxs'))
        self.assertEquals(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').split(), (u'/SNS', u'REF_M_1.nxs+REF_M_3.nxs'))

    def test_run_numbers(self):
        self.assertEquals(FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs').run_numbers(), [1, 3])
        file_path = FilePath(u'/SNS/REF_M_3.nxs+/SNS/REF_M_1.nxs+/SNS/REF_M_6.nxs+/SNS/REF_M_2.nxs')
        self.assertEqual(file_path.run_numbers(string_representation='long'), '1+2+3+6')
        self.assertEqual(file_path.run_numbers(string_representation='short'), '1:3+6')


if __name__ == '__main__':
    unittest.main()

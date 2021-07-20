from __future__ import absolute_import, division, print_function, unicode_literals

# package imports
from reflectivity_ui.interfaces.main_window import MainWindow
from reflectivity_ui.interfaces.event_handlers.main_handler import MainHandler

# 3rd-party imports
from PyQt5.QtWidgets import QApplication
import pytest

# standard imports
import os
import sys

this_module_path = sys.modules[__name__].__file__


class DataManagerMock(object):
    current_directory = os.path.dirname(this_module_path)


class MainWindowMock(object):
    ui = None
    main_window = None
    data_manager = DataManagerMock()


class TestMainHandler(object):

    app = QApplication(sys.argv)
    application = MainWindow()
    handler = MainHandler(application)

    def test_congruency_fail_report(self, data_server):
        # Selected subset of log names with an invalid one
        message = self.handler._congruency_fail_report([data_server.path_to('REF_M_24945_event.nxs'),
                                                        data_server.path_to('REF_M_24949_event.nxs')],
                                                       log_names=['LambdaRequest', 'NoLog'])
        assert message == 'NoLog is not a valid Log for comparison'

        # Valid subset of log name
        message = self.handler._congruency_fail_report([data_server.path_to('REF_M_24945_event.nxs'),
                                                        data_server.path_to('REF_M_24949_event.nxs')],
                                                       log_names=['LambdaRequest', 'frequency'])
        assert message == ''

        # Old files
        message = self.handler._congruency_fail_report([data_server.path_to('REF_M_24945_event.nxs'),
                                                        data_server.path_to('REF_M_24949_event.nxs')])
        assert 'values for log S3Vheight that differ above tolerance 0.01' in message

        # New files
        message = self.handler._congruency_fail_report([data_server.path_to('REF_M_38198.nxs.h5'),
                                                        data_server.path_to('REF_M_38199.nxs.h5')])
        assert 'values for log S1HWidth that differ above tolerance 0.01' in message


if __name__ == '__main__':
    pytest.main([__file__])

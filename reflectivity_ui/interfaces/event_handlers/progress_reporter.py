"""
    Class used to report on progress. It allows for sub-tasks and
    computes a meaningful progress status accordingly.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
from PyQt5 import QtGui, QtCore, QtWidgets

class ProgressReporter(object):

    def __init__(self, max_value=100, call_back=None, create_dialog=False,
                 window_title="Loading", message="", parent=None):
        self.max_value = max_value
        self.call_back = call_back
        self.value = 0
        self.sub_tasks = []
        self.main_window = parent
        self.progress_dialog = None

        if create_dialog:
            self.progress_dialog = QtWidgets.QProgressDialog(
                message, "Close", 0, max_value, self.main_window)
            self.progress_dialog.setWindowTitle(window_title)
            self.progress_dialog.setAutoClose(True)
            self.progress_dialog.setModal(True)
            self.progress_dialog.show()
            QtWidgets.QApplication.instance().processEvents()

    def __call__(self, value, message='', out_of=None):
        """
            Shortcut to set_value() so that the object can be used
            as a function to be compatible with QProgressDialog.setValue().
        """
        return self.set_value(value, message, out_of)

    def set_value(self, value, message='', out_of=None):
        """
            Set the value of a progress indicator
            :param int value: completion value, as a percentage
        """
        if out_of is not None:
            value = int(value / out_of * self.max_value)
        value = min(value, self.max_value)
        self.value = value
        self.update(message)

    def update(self, message=''):
        """
            Updates the progress status according to
            sub-tasks.
        """
        _value = self.value
        for item in self.sub_tasks:
            _value += min(item.value, item.max_value)
        _value = min(_value, self.max_value)

        if self.call_back is not None:
            self.call_back(message)
        if self.progress_dialog:
            self.progress_dialog.setValue(_value)
            self.progress_dialog.setLabelText(message)
            QtWidgets.QApplication.instance().processEvents()

    def create_sub_task(self, max_value):
        """
            Create a sub-task, with max_value being its portion
            of the complete task. Returns a call-back function
            to be called by the worker to update the progress.

            :param int max_value: portion of the task
        """
        p = ProgressReporter(max_value, self.update)
        self.sub_tasks.append(p)
        return p

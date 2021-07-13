"""
Zoo for customized simple widgets
"""
from PyQt5 import QDialog, QDialogButtonBox, QVBoxLayout, QLabel


class CustomDialog(QDialog):
    """
    Refer to: https://www.mfitzp.com/tutorials/pyqt-dialogs/
    """
    def __init__(self, parent=None, title='', message=''):
        QDialog.__init__(parent)

        self.setWindowTitle("HELLO!")

        QBtn = QDialogButtonBox.Ok | QDialogButtonBox.Cancel

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)

        self.layout = QVBoxLayout()
        message = QLabel("Something happened, is that OK?")
        self.layout.addWidget(message)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

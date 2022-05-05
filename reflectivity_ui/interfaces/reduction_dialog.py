"""
   Dialog to select reduction options to choose which outputs are needed
   and in which formats to write them.
"""
#pylint: disable=bare-except

import os
from PyQt5 import QtCore, QtWidgets
import reflectivity_ui.interfaces.generated.ui_reduce_dialog

class ReductionDialog(QtWidgets.QDialog, reflectivity_ui.interfaces.generated.ui_reduce_dialog.Ui_Dialog):
    """
        Reduction dialog
    """
    default_template = '(instrument)_{numbers}_{item}_{state}.{type}'

    def __init__(self, parent):
        super(ReductionDialog, self).__init__(parent)

        self.setupUi(self)

        self.settings = QtCore.QSettings('.refredm')

        self.directoryEntry.setText(self.settings.value('output_directory',
                                                        os.path.expanduser('~')))
        self.fileNameEntry.setText(self.settings.value('output_file_template',
                                                       ReductionDialog.default_template))

        # Outputs
        self.exportSpecular.setChecked(self._verify_true('export_specular', True))
        self.export_SA.setChecked(self._verify_true('export_asym', True))
        self.exportGISANS.setChecked(self._verify_true('export_gisans', False))
        self.exportOffSpecular.setChecked(self._verify_true('export_offspec', False))
        self.exportOffSpecularSmoothed.setChecked(self._verify_true('export_offspec_smooth', False))

        # Formats
        self.genx.setChecked(self._verify_true('format_genx', False))
        self.matlab.setChecked(self._verify_true('format_matlab', False))
        self.multiAscii.setChecked(self._verify_true('format_multi', False))
        self.numpy.setChecked(self._verify_true('format_numpy', False))
        self.mantid_script_checkbox.setChecked(self._verify_true('format_mantid', False))
        self.five_cols_checkbox.setChecked(self._verify_true('format_5cols', True))

        # Emails
        self.emailSend.setChecked(self._verify_true('email_send', False))
        self.emailZIPData.setChecked(self._verify_true('email_zip_data', True))
        self.emailSendPlots.setChecked(self._verify_true('email_send_plots', False))
        self.emailSendData.setChecked(self._verify_true('email_send_data', True))
        self.emailTo.setText(self.settings.value('email_to', ''))
        self.emailCc.setText(self.settings.value('email_cc', ''))
        self.emailSubject.setText(self.settings.value('email_subject', ''))

        self.is_accepted = False

    def _verify_true(self, parameter, default):
        """ Utility function to read a bool """
        _value = self.settings.value(parameter, str(default))
        return str(_value).lower() == 'true'

    def accept(self):
        """
            Save the current options and close dialog
        """
        self.save_settings()
        self.is_accepted = True
        self.close()

    def get_options(self):
        """ Return the reduction options as a dict"""
        if self.is_accepted is False:
            return None
        return dict(export_specular=self.exportSpecular.isChecked(),
                    export_asym=self.export_SA.isChecked(),
                    export_gisans=self.exportGISANS.isChecked(),
                    export_offspec=self.exportOffSpecular.isChecked(),
                    export_offspec_smooth=self.exportOffSpecularSmoothed.isChecked(),
                    format_genx=self.genx.isChecked(),
                    format_matlab=self.matlab.isChecked(),
                    format_mantid=self.mantid_script_checkbox.isChecked(),
                    format_multi=self.multiAscii.isChecked(),
                    format_numpy=self.numpy.isChecked(),
                    format_5cols=self.five_cols_checkbox.isChecked(),
                    output_directory=self.directoryEntry.text(),
                    output_file_template=self.fileNameEntry.text(),
                    email_send=self.emailSend.isChecked(),
                    email_zip_data=self.emailZIPData.isChecked(),
                    email_send_plots=self.emailSendPlots.isChecked(),
                    email_send_data=self.emailSendData.isChecked(),
                    email_to=self.emailTo.text(),
                    email_cc=self.emailCc.text(),
                    email_subject=self.emailSubject.text(),
                    email_text=self.emailText.toPlainText(),
                    )

    def change_directory(self):
        """
            Change the output directory
        """
        old_d = self.directoryEntry.text()
        new_d = QtWidgets.QFileDialog.getExistingDirectory(parent=self,
                                                           caption='Select new directory',
                                                           directory=old_d)
        if new_d is not None:
            self.directoryEntry.setText(new_d)

    def save_settings(self):
        """
            Save reduction options in QSettings
        """
        self.settings.setValue('output_directory', self.directoryEntry.text())
        self.settings.setValue('output_file_template', self.fileNameEntry.text())

        self.settings.setValue('export_specular', self.exportSpecular.isChecked())
        self.settings.setValue('export_asym', self.export_SA.isChecked())
        self.settings.setValue('export_gisans', self.exportGISANS.isChecked())
        self.settings.setValue('export_offspec', self.exportOffSpecular.isChecked())
        self.settings.setValue('export_offspec_smooth', self.exportOffSpecularSmoothed.isChecked())

        self.settings.setValue('format_genx', self.genx.isChecked())
        self.settings.setValue('format_matlab', self.matlab.isChecked())
        self.settings.setValue('format_multi', self.multiAscii.isChecked())
        self.settings.setValue('format_numpy', self.numpy.isChecked())
        self.settings.setValue('format_mantid', self.mantid_script_checkbox.isChecked())
        self.settings.setValue('format_5cols', self.five_cols_checkbox.isChecked())

        self.settings.setValue('email_send', self.emailSend.isChecked())
        self.settings.setValue('email_zip_data', self.emailZIPData.isChecked())
        self.settings.setValue('email_send_plots', self.emailSendPlots.isChecked())
        self.settings.setValue('email_send_data', self.emailSendData.isChecked())
        self.settings.setValue('email_to', self.emailTo.text())
        self.settings.setValue('email_cc', self.emailCc.text())
        self.settings.setValue('email_subject', self.emailSubject.text())

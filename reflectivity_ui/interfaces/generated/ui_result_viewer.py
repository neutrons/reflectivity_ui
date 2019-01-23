# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer/ui_result_viewer.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(1031, 552)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget = QtWidgets.QTabWidget(Dialog)
        self.tabWidget.setObjectName("tabWidget")
        self.specular_tab = QtWidgets.QWidget()
        self.specular_tab.setObjectName("specular_tab")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.specular_tab)
        self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.specular_compare_widget = CompareWidget(self.specular_tab)
        self.specular_compare_widget.setObjectName("specular_compare_widget")
        self.verticalLayout_2.addWidget(self.specular_compare_widget)
        self.tabWidget.addTab(self.specular_tab, "")
        self.offspecular_tab = QtWidgets.QWidget()
        self.offspecular_tab.setObjectName("offspecular_tab")
        self.tabWidget.addTab(self.offspecular_tab, "")
        self.gisans_tab = QtWidgets.QWidget()
        self.gisans_tab.setObjectName("gisans_tab")
        self.tabWidget.addTab(self.gisans_tab, "")
        self.verticalLayout.addWidget(self.tabWidget)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        self.tabWidget.setCurrentIndex(0)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.specular_tab), _translate("Dialog", "Specular"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.offspecular_tab), _translate("Dialog", "Off-Specular"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.gisans_tab), _translate("Dialog", "GISANS"))

from .compare_plots import CompareWidget

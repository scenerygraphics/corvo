from functools import partial
from os import listdir
from os.path import isfile, join

from PyQt5.QtCore import pyqtSlot, Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QLabel, QPushButton, QMessageBox

from corvolauncher.gui.qt_line_break import QHLineBreakWidget
from corvolauncher.gui.qt_process_menu import ProcessMenu
from corvolauncher.utilities.jar_launch import JarLaunch


class LaunchWindow(QWidget):
    def __init__(self, parent, dataset):
        super().__init__(parent)

        self.parent = parent
        self.dataset = dataset

        msg = QMessageBox()
        msg.setWindowTitle("Corvo Launcher")
        msg.setText("Launch Corvo with " + self.dataset + "?")
        msg.setIcon(QMessageBox.Question)
        msg.setStandardButtons(QMessageBox.Cancel | QMessageBox.Yes)
        msg.setDefaultButton(QMessageBox.Yes)

        msg.buttonClicked.connect(self.launch)
        msg.exec()

    def launch(self, button):
        if button.text() == "&Yes":
            print(self.dataset)

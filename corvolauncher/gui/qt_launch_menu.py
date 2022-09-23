import threading

from PyQt5.QtCore import pyqtSlot
from PyQt5.QtWidgets import QHBoxLayout, QLabel, QPushButton, QWidget

from corvolauncher.utilities.jar_launch import JarLaunch


class LaunchMenu(QWidget):
    def __init__(self, parent, dataset):
        super().__init__()
        self.parent = parent
        self.dataset = dataset

        self.layout = QHBoxLayout()
        self.setFixedWidth(2 * self.parent.win_width // 3)

        self.launch_button = QPushButton("Launch " + self.dataset)
        self.launch_button.clicked.connect(self.launch_corvo)
        self.layout.addWidget(self.launch_button)
        self.setLayout(self.layout)
        # self.show()

    @pyqtSlot()
    def launch_corvo(self):
        # print("launch")
        # JarLaunch(self.dataset)
        pass
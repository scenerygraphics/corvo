from functools import partial
from os import listdir
from os.path import isfile, join

from PyQt5.QtCore import pyqtSlot, Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QLabel, QPushButton, QMessageBox

from corvolauncher.gui.qt_line_break import QHLineBreakWidget
from corvolauncher.gui.qt_process_menu import ProcessMenu
from corvolauncher.gui.qt_launch_window import LaunchWindow


class DatasetSidebar(QWidget):
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool
        self.current_popup = ""

        self.title_font = QFont()
        self.title_font.setBold(True)

        self.master_layout = QGridLayout()
        self.master_layout.setVerticalSpacing(20)
        self.master_layout.setAlignment(Qt.AlignTop)
        self.setLayout(self.master_layout)
        self.setFixedWidth(self.parent.win_width // 3)

        self.raw_layout = QVBoxLayout()
        self.processed_layout = QVBoxLayout()

        self.raw_layout.setSpacing(5)
        self.processed_layout.setSpacing(5)

        self.raw_title = QLabel("Raw Datasets")
        self.processed_title = QLabel("Processed Datasets")

        for t in [self.raw_title, self.processed_title]:
            t.setFont(self.title_font)

        self.raw_layout.addWidget(self.raw_title)
        self.processed_layout.addWidget(self.processed_title)

        self.raw_datasets = [f for f in listdir("../resources/datasets/") if isfile(join("../resources/datasets/", f))]

        for file_name in self.raw_datasets:
            # partial allows connection to be made to declaration, not overriding connecting object as with lambda
            button = QPushButton(file_name)
            self.raw_layout.addWidget(button)
            button.clicked.connect(partial(self.launch_config, button.text()))

        self.processed_datasets = [f for f in listdir("../resources/processed_datasets/") if
                                   isfile(join("../resources/processed_datasets/", f))]

        for file_name in self.processed_datasets:
            button = QPushButton(file_name)
            self.processed_layout.addWidget(button)
            button.clicked.connect(partial(self.launch_popup, button.text()))

        self.master_layout.addLayout(self.raw_layout, 0, 0)
        self.master_layout.addWidget(QHLineBreakWidget(self), 1, 0)
        self.master_layout.addLayout(self.processed_layout, 2, 0)
        self.master_layout.addWidget(QHLineBreakWidget(self), 3, 0)

    @pyqtSlot(str)
    def launch_config(self, button):
        self.parent.file_tabs.addTab(ProcessMenu(self, button, self.threadpool), button)

    @pyqtSlot(str)
    def launch_popup(self, dataset):
        LaunchWindow(self, dataset)

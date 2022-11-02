import os
from functools import partial
from os import listdir
from os.path import isfile, join

from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QWidget, QGridLayout, QVBoxLayout, QLabel, QPushButton, QHBoxLayout, QFrame, \
    QGraphicsDropShadowEffect, QSizePolicy

from corvolauncher.gui.dataset_fetcher import DatasetFetcher
from corvolauncher.gui.job_runners.generic_worker import GenericWorker
from corvolauncher.gui.job_runners.jar_worker import JarWorker
from corvolauncher.gui.qt_confirm_popup import ConfirmPopup
from corvolauncher.gui.qt_line_break import QHLineBreakWidget
from corvolauncher.gui.qt_process_menu import ProcessMenu


class DatasetSidebar(QWidget):
    trigger_stop_corvo = pyqtSignal()
    # corvo_running_flag = False
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        #  create and configure sidebar master and child layob nuts
        # self.master_layout = QGridLayout()
        self.master_layout = QVBoxLayout()
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # self.master_layout.setVerticalSpacing(10)
        self.master_layout.setSpacing(10)
        self.master_layout.setAlignment(Qt.AlignTop)


        self.setLayout(self.master_layout)
        self.setFixedWidth(self.parent.win_width * 3 // 4)
        # self.setFixedWidth(self.parent.win_width // 3)

        self.raw_layout = QVBoxLayout()
        self.processed_layout = QVBoxLayout()
        self.raw_layout.setSpacing(5)
        self.processed_layout.setSpacing(5)

        #  give section titles a bold font
        self.raw_title = QLabel("Raw Datasets")
        self.processed_title = QLabel("Processed Datasets")
        self.title_font = QFont()
        self.title_font.setPointSize(10)
        self.title_font.setBold(True)
        for t in [self.raw_title, self.processed_title]:
            t.setFont(self.title_font)

        self.raw_layout.addWidget(self.raw_title)
        self.processed_layout.addWidget(self.processed_title)

        # print(self.processed_layout.count())
        # print(self.processed_layout.itemAt(0).widget().objectName())

        #  create shutdown button for new Corvo instance
        self.close_corvo_button = QPushButton("shut down Corvo")
        self.close_corvo_button.hide()
        self.close_corvo_button.clicked.connect(self.shutdown_corvo)
        self.processed_layout.addWidget(self.close_corvo_button)

        self.add_raw_files([f for f in listdir("../resources/datasets/") if isfile(join("../resources/datasets/", f))])
        self.add_processed_files([f for f in listdir("../resources/processed_datasets/") if
                                  isfile(join("../resources/processed_datasets/", f))])

        self.master_layout.addLayout(self.raw_layout)
        self.master_layout.addWidget(QHLineBreakWidget(self))
        self.master_layout.addLayout(self.processed_layout)
        self.master_layout.addWidget(QHLineBreakWidget(self))
        self.dataset_fetcher = DatasetFetcher(self, self.threadpool)
        self.master_layout.addWidget(self.dataset_fetcher)
        self.dataset_fetcher.collection_container.blockSignals(False)
        self.master_layout.addStretch()

        # self.master_layout.addLayout(self.raw_layout, 0, 0)
        # self.master_layout.addWidget(QHLineBreakWidget(self), 1, 0)
        # self.master_layout.addLayout(self.processed_layout, 2, 0)
        # self.master_layout.addWidget(QHLineBreakWidget(self), 3, 0)
        # self.dataset_fetcher = DatasetFetcher(self, self.threadpool)
        # self.master_layout.addWidget(self.dataset_fetcher)
        # self.dataset_fetcher.collection_container.blockSignals(False)


    @pyqtSlot(str)
    def rm_file(self, file_name):
        @pyqtSlot()
        def on_finished():
            # slot triggered on finish to update sidebar
            self.update_directory_index()

        def rm_os(f_name):
            #  function to be encapsulated in pop-up window
            def container():
                try:
                    os.remove("../resources/datasets/" + f_name)
                except FileNotFoundError:
                    os.remove("../resources/processed_datasets/" + f_name)

            worker = GenericWorker(container)
            worker.signals.finished.connect(on_finished)
            self.threadpool.start(worker)

        ConfirmPopup("Are you sure you want to delete this dataset from your system?", rm_os, file_name)

    @pyqtSlot()
    def shutdown_corvo(self):
        self.trigger_stop_corvo.emit()

    @pyqtSlot(str)
    def launch_jar(self, file_name):
        @pyqtSlot()
        def on_running():
            #  disable dataset launch buttons to prevent duplicate launches
            self.processed_layout.itemAt(2).widget().setEnabled(False)
            self.processed_layout.itemAt(1).widget().show()

        @pyqtSlot()
        def on_finished():
            self.processed_layout.itemAt(2).widget().setEnabled(True)
            self.processed_layout.itemAt(1).widget().hide()

        def container(f_name):
            worker = JarWorker(f_name)
            worker.signals.running.connect(on_running)
            worker.signals.finished.connect(on_finished)
            self.trigger_stop_corvo.connect(worker.stop)

            self.threadpool.start(worker)

        ConfirmPopup("Launch Corvo with " + file_name[:-31].replace("_", " ") + "?", container, file_name)

    def update_directory_index(self):
        self.raw_layout.removeWidget(self.raw_layout.itemAt(1).widget())  # removing item doesnt seem to work - parse layout parent instead
        self.processed_layout.removeWidget(self.processed_layout.itemAt(2).widget())

        self.add_raw_files([f for f in listdir("../resources/datasets/") if isfile(join("../resources/datasets/", f))])
        self.add_processed_files([f for f in listdir("../resources/processed_datasets/") if
                                isfile(join("../resources/processed_datasets/", f))])

    def add_raw_files(self, r_files: list):
        layout_container = QFrame()
        button_layout = QVBoxLayout()
        for file_name in r_files:
            if file_name != "":  # potentially a .DS_Store file on mac
                # partial allows connection to be made to declaration, not overriding connecting object as with lambda
                button_del_set = QHBoxLayout()
                button_del_set.setAlignment(Qt.AlignLeft)

                button_del_set.setEnabled(False)

                button = QPushButton("pre-process")
                button.setStyleSheet(":enabled{background-color:rgb(201, 132, 46);color: black;font-size: 9pt;}" "QPushButton:pressed {background-color:rgb(82, 64, 33);color: black;font-size: 9pt;}")
                shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=3, yOffset=3)
                button.setGraphicsEffect(shadow)
                button.setObjectName(file_name)
                button.clicked.connect(partial(self.configuration_window, button.objectName()))

                del_button = QPushButton("x")
                del_button.setStyleSheet("background-color:rgb(145, 45, 45)")
                del_button.clicked.connect(partial(self.rm_file, button.objectName()))

                dset_label = QLabel(file_name[:-15].replace("_", " "))
                dset_label.setFixedWidth(self.width() - 125)
                dset_label.setWordWrap(True)

                button_del_set.addWidget(del_button)
                button_del_set.addWidget(button)
                button_del_set.addWidget(dset_label)
                button_layout.addLayout(button_del_set)

        layout_container.setLayout(button_layout)
        layout_container.setObjectName("raw_frame")
        self.raw_layout.addWidget(layout_container)

    def add_processed_files(self, p_files: list):
        layout_container = QFrame()
        button_layout = QVBoxLayout()
        for file_name in p_files:
            if file_name != "":
                button_del_pair = QHBoxLayout()
                button_del_pair.setAlignment(Qt.AlignLeft)

                button_del_pair.setEnabled(False)

                button = QPushButton("launch")
                button.setStyleSheet(":enabled{background-color:rgb(67, 171, 112);color: black;font-size: 9pt;}" "QPushButton:pressed {background-color:rgb(36, 99, 63);color: black;font-size: 9pt;}")
                shadow = QGraphicsDropShadowEffect(blurRadius=5, xOffset=3, yOffset=3)
                button.setGraphicsEffect(shadow)
                button.setObjectName(file_name)
                button.clicked.connect(partial(self.launch_jar, button.objectName()))

                del_button = QPushButton("x")
                del_button.setStyleSheet("background-color:rgb(145, 45, 45)")
                del_button.clicked.connect(partial(self.rm_file, button.objectName()))

                dset_label = QLabel(file_name[:-31].replace("_", " "))
                dset_label.setFixedWidth(self.width() - 90)
                dset_label.setWordWrap(True)

                button_del_pair.addWidget(del_button)
                button_del_pair.addWidget(button)
                button_del_pair.addWidget(dset_label)
                button_layout.addLayout(button_del_pair)

        layout_container.setObjectName("processed_frame")
        layout_container.setLayout(button_layout)
        self.processed_layout.addWidget(layout_container)

    @pyqtSlot(str)
    def configuration_window(self, button):
        # make sure a tab with this dataset is not currently open
        self.parent.file_tabs.addTab(ProcessMenu(self, button, self.threadpool), button[:-14].replace("_", " ")[0:20] + "...")

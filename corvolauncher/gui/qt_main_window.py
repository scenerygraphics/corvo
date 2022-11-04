import json
import os.path
import sys
import threading
import time

import qdarkstyle
from pathlib import Path

# from PyQt5.QtCore import Qt
from PyQt5.QtCore import Qt, QThreadPool, pyqtSlot
from PyQt5.QtWidgets import (
    QMainWindow,
    QStatusBar,
    QApplication,
    QWidget,
    QDockWidget, QGridLayout, QHBoxLayout, QLabel, QTabWidget, QTabBar,
)

from corvolauncher.gui.qt_sidebar import DatasetSidebar


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title = "Corvo Launcher"
        self.threadpool = QThreadPool()

        self.desktop = QApplication.desktop()
        self.screenRect = self.desktop.screenGeometry()
        height, width = self.screenRect.height(), self.screenRect.width()
        self.setAcceptDrops(True)

        self.win_width = width // 3
        self.height = height // 2
        # self.height = 650
        self.left = (width - self.win_width) // 2
        self.top = (height - self.height) // 2

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.win_width, self.height)

        self.setDockOptions(QMainWindow.AnimatedDocks | QMainWindow.AllowNestedDocks)
        # self.setDockOptions(QMainWindow.ForceTabbedDocks)

        self.file_display_dock = QDockWidget(self)

        self.file_display_dock.setWindowTitle("Compatible Files")
        self.file_display_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)

        self.side_bar = DatasetSidebar(self, self.threadpool)

        self.file_display_dock.setWidget(self.side_bar)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.file_display_dock)

        self.file_tabs = QTabWidget(self, movable=True, tabsClosable=True)
        self.file_tabs.tabCloseRequested.connect(self.close_tab)
        self.no_dataset_label = QLabel("Drag and drop your .h5ad dataset into the window")
        self.no_dataset_label.setAlignment(Qt.AlignCenter)
        self.file_tabs.addTab(self.no_dataset_label, "")
        self.file_tabs.tabBar().setTabButton(0, QTabBar.RightSide, None)

        self.setCentralWidget(self.file_tabs)

        self.show()

    # def add_as_dock(self, widget: QWidget, title: str):
    #     dock = QDockWidget(self)
    #     dock.setWindowTitle(title + " setup")
    #     dock.setFeatures(QDockWidget.AllDockWidgetFeatures)
    #     dock.setAllowedAreas(Qt.RightDockWidgetArea)
    #     dock.setWidget(widget)
    #     self.addDockWidget(Qt.RightDockWidgetArea, dock)

    @pyqtSlot(int)
    def close_tab(self, index):
        self.file_tabs.removeTab(index)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        for f in files:
            d_name = f.split("/")[-1]
            d_type = f[-5:]
            if d_type == ".h5ad" and not f[-20:-5] == "corvo_PROCESSED" and (
                    f[-14:-5] == "corvo_RAW" or f[-14:-7] == "corvo_RAW"):  # unprocessed h5ad datasets
                os.rename(f, "../resources/{}/{}".format("datasets", self.moniker_recursively(d_name, "datasets")))
            elif d_type == ".h5ad" and (
                    f[-20:-5] == "corvo_PROCESSED" or f[-20:-7] == "corvo_PROCESSED"):  # processed h5ad datasets
                os.rename(f, "../resources/{}/{}".format("processed_datasets",
                                                         self.moniker_recursively(d_name, "processed_datasets")))
            elif d_type == ".h5ad" and not f[-20:-5] == "corvo_PROCESSED" and not f[-14:-5] == "corvo_RAW":
                os.rename(f, "../resources/{}/{}".format("datasets", self.moniker_recursively(
                    d_name.rstrip(".h5ad") + "_corvo_RAW.h5ad", "datasets")))
            else:
                print("this file type is not supported")
        self.side_bar.update_directory_index()

    def moniker_recursively(self, file_name: str, dset_dir):
        if not os.path.exists("../resources/{}/{}".format(dset_dir, file_name)):
            return file_name
        else:
            num = file_name[-6]
            print(num)
            if num == "D" or num == "W":
                file_name = file_name.rstrip(".h5ad") + "_2" + ".h5ad"
                unique_name = self.moniker_recursively(file_name, dset_dir)
                return unique_name
            else:
                file_name = file_name.rstrip(num + ".h5ad") + str(int(num) + 1) + ".h5ad"
                unique_name = self.moniker_recursively(file_name, dset_dir)
                return unique_name


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    sys.exit(app.exec())

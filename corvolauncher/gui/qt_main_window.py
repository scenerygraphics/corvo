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

from corvolauncher.gui.qt_dataset_sidebar import DatasetSidebar


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
        # add recursive function to append an integer to the name of a file if the name already exists
        files = [u.toLocalFile() for u in event.mimeData().urls()]
        for f in files:
            d_name = f.split("/")[-1]
            d_type = f[-5:]
            if d_type == ".h5ad" and not os.path.exists("../resources/datasets/" + d_name) and not f[-14:-5] == "processed":
                os.rename(f, "../resources/datasets/" + d_name)
            elif d_type == ".h5ad" and not os.path.exists("../resources/processed_datasets/" + d_name) and f[-14:-5] == "processed":
                os.rename(f, "../resources/processed_datasets/" + d_name)
            else:
                print("the file type is not supported, or a dataset with the same name already exists")
        self.side_bar.update_directory_index()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    sys.exit(app.exec())

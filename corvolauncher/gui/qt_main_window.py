import json
import os.path
import sys
import threading
import time

import qdarkstyle
from pathlib import Path

# from PyQt5.QtCore import Qt
from PyQt5.QtCore import Qt, QThreadPool
from PyQt5.QtWidgets import (
    QMainWindow,
    QStatusBar,
    QApplication,
    QWidget,
    QDockWidget, QGridLayout, QHBoxLayout, QLabel,
)

from corvolauncher.gui.qt_dataset_sidebar import DatasetSidebar


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title = "Corvo Launcher"

        self.desktop = QApplication.desktop()
        self.screenRect = self.desktop.screenGeometry()
        height, width = self.screenRect.height(), self.screenRect.width()

        self.win_width = width // 3
        self.height = height // 2
        self.left = (width - self.win_width) // 2
        self.top = (height - self.height) // 2

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.win_width, self.height)

        # self.setDockOptions(QMainWindow.AnimatedDocks | QMainWindow.AllowNestedDocks)
        self.setDockOptions(QMainWindow.ForceTabbedDocks)

        self.file_display_dock = QDockWidget(self)

        self.file_display_dock.setWindowTitle("Compatible Files")
        self.file_display_dock.setFeatures(QDockWidget.NoDockWidgetFeatures)

        self.file_display_dock.setWidget(DatasetSidebar(self))
        self.addDockWidget(Qt.LeftDockWidgetArea, self.file_display_dock)

        self.setCentralWidget(QWidget())

        self.show()

    def add_as_dock(self, widget: QWidget, title: str):
        dock = QDockWidget(self)

        dock.setWindowTitle(title + " setup")

        dock.setFeatures(QDockWidget.AllDockWidgetFeatures)
        dock.setAllowedAreas(Qt.RightDockWidgetArea)

        dock.setWidget(widget)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    sys.exit(app.exec())

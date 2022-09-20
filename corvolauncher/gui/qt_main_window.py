import json
import os.path
import sys
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

from corvolauncher.gui.qt_file_display import FileDisplay


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

        self.layout = QHBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)
        # self.grid_layout = QGridLayout()
        # self.grid_layout.setAlignment(Qt.AlignTop)

        self.qt_file_display = FileDisplay(self)

        # self.grid_layout.addWidget(self.qt_file_display, 0, 0)
        self.layout.addWidget(self.qt_file_display)

        self.layout.addWidget(QLabel("test"))

        central_widget = QWidget()
        # central_widget.setLayout(self.grid_layout)
        central_widget.setLayout(self.layout)
        self.setCentralWidget(central_widget)

        self.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    sys.exit(app.exec())

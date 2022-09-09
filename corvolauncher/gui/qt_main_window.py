import json
import os.path
import sys
import qdarkstyle
from pathlib import Path

# from PyQt5.QtCore import Qt, QThreadPool
from PyQt5.QtWidgets import (
    QMainWindow,
    QStatusBar,
    QApplication,
    QWidget,
    QDockWidget,
)

from qt_file_display import FileDisplay


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.title = "Corvo Launcher"
        self.left = 10
        self.top = 10
        self.width = 900
        self.height = 1000

        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        # self.test_widget = QWidget(self)

        self.file_display_widget = FileDisplay(self)
        self.t_wid = QWidget()
        self.t_wid.setWindowTitle("test")

        self.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt5())
    sys.exit(app.exec())

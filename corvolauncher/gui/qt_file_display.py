from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QWidget, QGridLayout, QLabel


class FileDisplay(QWidget):
    def __init__(self, parent):
        super(QWidget, self).__init__(parent)

        self.parent = parent

        self.grid_layout = QGridLayout()
        self.setLayout(self.grid_layout)
        self.grid_layout.setAlignment(Qt.AlignTop)

        self.my_font = QFont()
        self.my_font.setBold(True)

        self.raw_datasets_label = QLabel("Raw Datasets")
        self.raw_datasets_label.setFont(self.my_font)

        self.raw_datasets_label.setText("test")

        self.grid_layout.addWidget(self.raw_datasets_label, 0, 0)

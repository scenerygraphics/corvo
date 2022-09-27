import pathlib
import scanpy as sc
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QCheckBox, QPushButton


class ProcessMenu(QWidget):
    # change to permanent QTabWidget, starting with message, "select a file to view processing options"
    # create in beginning in MainWindow, with dataset buttons connected to method to call method within ProcessMenu
    # in DatasetSidebar: self.parent.process_menu.add_tab(parent, dataset)
    # now with fixed sidebar, tabbed config windows, and popup launcher, move away from docks

    def __init__(self, parent, dataset: str):
        super().__init__()

        self.parent = parent
        self.dataset = dataset
        self.adata = sc.read_h5ad("../resources/datasets/" + self.dataset)

        self.layout = QVBoxLayout()
        self.layout.setAlignment(Qt.AlignTop)
        self.setFixedWidth(self.parent.width() * 2)
        self.setLayout(self.layout)
        #
        # for observation in self.adata.obs:
        #     self.layout.addWidget(QCheckBox(observation))

        self.process_button = QPushButton("Pre-process " + self.dataset)
        self.process_button.clicked.connect(self.handle_launch)
        self.layout.addWidget(self.process_button)

    def handle_launch(self):
        print("launching")
        pass

    def launch(self):
        pass

if __name__ == "__main__":
    ProcessMenu("marrow.h5ad")

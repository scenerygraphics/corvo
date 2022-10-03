import pathlib
import time

import scanpy as sc
import numpy as np
from PyQt5.QtCore import Qt, pyqtSlot
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QCheckBox, QPushButton

from corvolauncher.gui.worker import Worker
from corvolauncher.utilities.complete_process import CompleteProcess


class ProcessMenu(QWidget):
    # change to permanent QTabWidget, starting with message, "select a file to view processing options"
    # create in beginning in MainWindow, with dataset buttons connected to method to call method within ProcessMenu
    # in DatasetSidebar: self.parent.process_menu.add_tab(parent, dataset)
    # now with fixed sidebar, tabbed config windows, and popup launcher, move away from docks

    def __init__(self, parent, dataset: str, threadpool):
        super().__init__()

        self.stateTracker = False

        self.parent = parent
        self.dataset = dataset
        self.threadpool = threadpool
        # self.adata = sc.read_h5ad("../resources/datasets/" + self.dataset)

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
        print("handling launch of " + self.dataset)
        self.process_button.setDisabled(True)
        self.stateTracker = not self.stateTracker  # implement color scheme hints and process cancel
        self.launch()

    def launch(self):
        print("launching " + self.dataset)
        # process_worker = Worker(self.scanpyProcess(self.dataset))
        process_worker = Worker(self.wait_func())
        process_worker.signals.finished.connect(self.resetButton)
        self.threadpool.start(process_worker)

    @pyqtSlot()
    def resetButton(self):
        self.process_button.setEnabled(True)

    def wait_func(self):
        time.sleep(100)

    def scanpyProcess(self, dataset):
        print("scanpy processing " + dataset)
        self.dataset: str = dataset

        adata = sc.read_h5ad("../resources/datasets/" + dataset)
        # adata.var_names = adata.var["feature_name"]

        # file to add 3d umap to
        out_file: str = dataset.rstrip(".h5ad") + "_processed.h5ad"

        # will fail for obs with only one category or when a category has only one data point
        n_genes = 10

        for observation in adata.obs:
            try:

                sc.tl.rank_genes_groups(adata, groupby=observation, n_genes=n_genes, method='t-test', use_raw=True)

                print("success: " + observation)
                names_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                pvals_list = [0] * adata.obs[observation].cat.categories.size * n_genes
                logfoldchanges_list = [0] * adata.obs[observation].cat.categories.size * n_genes

                for cat in [["names", names_list], ["pvals", pvals_list], ["logfoldchanges", logfoldchanges_list]]:

                    rank_counter = 0
                    for rank in adata.uns["rank_genes_groups"][cat[0]]:  # rank is list of 1st, 2nd, 3d most expr etc
                        name_counter = 0
                        if cat[0] == "pvals":
                            for index in rank:
                                cat[1][(name_counter * n_genes) + rank_counter] = round(index, 4)
                                name_counter += 1

                        else:
                            for index in rank:
                                cat[1][(name_counter * n_genes) + rank_counter] = index
                                name_counter += 1
                                # print(index)

                        rank_counter += 1

                    adata.uns[(observation + "_" + cat[0])] = cat[1]

            # pass for any obs that do not cluster
            except (AttributeError, ValueError, ZeroDivisionError):
                print("fail: " + observation)
                pass

        # adata.obsm['X_umap_2d'] = adata.obsm['X_umap'].copy()

        # may need flattening
        try:
            adata.uns['neighbors']['params']['n_pcs'] = adata.uns['neighbors']['params']['n_pcs'][0]
        except KeyError:
            pass

        sc.tl.umap(adata, n_components=3)

        # making X sparse in csc format
        adata.X = adata.X.tocsc()

        adata.write(pathlib.Path("../resources/processed_datasets/" + out_file))

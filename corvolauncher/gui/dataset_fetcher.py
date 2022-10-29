import sys
import threading
from functools import partial

from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal, QStringListModel, QSortFilterProxyModel
from PyQt5.QtGui import QFont, QMouseEvent, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QWidget, QComboBox, QMenuBar, QAction, QHBoxLayout, QPushButton, \
    QListView, QCompleter, QLineEdit, QFrame, QApplication
from requests import Response

from corvolauncher.gui.job_runners.cellxgene_worker import CellXGeneDownloadWorker, CellXGeneJSONWorker
from corvolauncher.gui.job_runners.generic_worker import GenericWorker
from corvolauncher.gui.qt_line_break import QHLineBreakWidget
from corvolauncher.utilities.cellxgene_JSON_requests import CellxGeneJSONRequests


class DatasetFetcher(QWidget):
    cancel_download = pyqtSignal()
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        self.currently_downloading = False
        self.currently_loading = False

        self.fetcher_layout = QVBoxLayout()
        self.setLayout(self.fetcher_layout)
        self.fetcher_layout.setSpacing(10)
        self.fetcher_layout.setAlignment(Qt.AlignLeft)

        #  give section titles a bold font
        self.title = QLabel("CZCellxGene Hosted Datasets")
        self.title_font = QFont()
        self.title_font.setPointSize(10)
        self.title_font.setBold(True)
        self.title.setFont(self.title_font)
        self.fetcher_layout.addWidget(self.title)

        self.data_fetcher = CellxGeneJSONRequests(self)
        self.collection_map = self.data_fetcher.collections
        self.collection_keys = self.collection_map.keys()
        self.index_list = list(self.collection_keys)
        self.index_list.insert(0, "-- search for a collection --")

        # block changed signal emitted on creation to stop addition to layout that doesn't exist yet
        self.combo_font = QFont()
        self.combo_font.setPointSize(10)

        self.collection_search = QLineEdit("-- search for a collection --")
        self.collection_completer = QCompleter(self.collection_keys, self.collection_search)
        self.collection_completer.setCaseSensitivity(False)
        self.collection_search.setCompleter(self.collection_completer)
        self.collection_search.setFont(self.combo_font)
        self.collection_search.setFixedHeight(25)

        self.collection_search.editingFinished.connect(self.search_box_handler)
        self.fetcher_layout.addWidget(self.collection_search)

        self.collection_container = QComboBox()
        self.collection_container.blockSignals(True)
        self.collection_container.currentIndexChanged.connect(self.load_collection_dsets)
        self.collection_container.setFont(self.combo_font)
        self.collection_container.setFixedHeight(25)

        listView = QListView()
        listView.setWordWrap(True)
        listView.setSpacing(5)
        self.collection_container.setView(listView)
        self.collection_container.addItem("-- select a collection to load its datasets --")

        for c in self.collection_keys:
            self.collection_container.addItem(c)
        self.fetcher_layout.addWidget(self.collection_container)

    @pyqtSlot()
    def search_box_handler(self):
        try:
            collection_index = self.index_list.index(self.collection_search.text())
            self.load_collection_dsets(collection_index)
        except Exception:
            print("collection does not exist")
        # self.load_collection_dsets(self.index_list.index(self.collection_search.text()))

    @pyqtSlot(int)
    def load_collection_dsets(self, index):
        # self.currently_loading = True
        print("trying to load index: " + str(index))
        # loading datasets still runs blocking and overruns the size of the window if too many datasets exist
        if self.fetcher_layout.count() > 3:
            self.fetcher_layout.removeWidget(self.children()[4])


        print(self.fetcher_layout.count())
        print(len(self.children()))
        if index != 0:
            dataset_map = self.data_fetcher.get_collection_datasets(
                self.collection_map[self.collection_container.itemText(index)])

            dset_button_layout = QVBoxLayout()
            dset_button_layout.setSpacing(10)
            dset_button_frame = QFrame()
            dset_button_layout.setAlignment(Qt.AlignLeft)

            for dset in dataset_map.keys():
                dset_layout = QHBoxLayout()
                dset_download_button = QPushButton("Download")
                dset_download_button.setObjectName(dset)
                dset_download_button.setFixedWidth(65)
                dset_download_button.setFixedHeight(20)
                dset_download_button.clicked.connect(partial(self.launch_download, dset, dataset_map))
                dset_layout.addWidget(dset_download_button)

                # dset_cancel_button = QPushButton("Cancel")
                # dset_cancel_button.setObjectName(dset + "cancel")
                # dset_download_button.setFixedWidth(65)
                # dset_cancel_button.clicked.connect(self.cancel)
                # dset_cancel_button.hide()
                # dset_layout.addWidget(dset_cancel_button)

                dset_label = QLabel(dset)
                dset_label.setMinimumHeight(30)
                dset_label.setFixedWidth(self.collection_container.width() - dset_download_button.width())
                # dset_label.setFixedWidth(self.collection_container.width() - 125)
                # dset_label.setMaximumWidth(self.collection_container.width() - dset_download_button.width())
                dset_label.setWordWrap(True)
                dset_layout.addWidget(dset_label)
                dset_button_layout.addLayout(dset_layout)
            dset_button_frame.setLayout(dset_button_layout)
            self.fetcher_layout.addWidget(dset_button_frame)

        # @pyqtSlot()
        # def on_running():
        #     self.fetcher_layout.addWidget(QLabel("fetching..."))
        #
        # @pyqtSlot()
        # def on_finished():
        #     self.fetcher_layout.removeWidget(self.children()[3])

        # worker = CellXGeneJSONWorker(self, index)
        # self.threadpool.start(worker)
        # worker = GenericWorker(container, index)
        # worker.signals.running.connect(on_running)
        # worker.signals.finished.connect(on_finished)
        # self.threadpool.start(worker)

    # @pyqtSlot()
    # def cancel(self):
    #     self.cancel_download.emit()

    @pyqtSlot(str, dict)
    def launch_download(self, name, ids):
        print(name)
        print(ids)

        # btn_pressed = self.sender()
        # btn_pressed.hide()
        # @pyqtSlot()
        # def on_running():
        #     btn_pressed.hide()
        #     print(self.fetcher_layout.findChild(QPushButton, name + "cancel"))

        @pyqtSlot()
        def on_finished():
            self.parent.update_directory_index()

        @pyqtSlot(float)
        def disp_file_size(size):
            print(float(size) / 1000000)

        @pyqtSlot(int)
        def progress_bar(done):
            sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50 - done)))
            sys.stdout.flush()

        # # worker = GenericWorker(self.data_fetcher.download_dataset, ids, name)
        worker = CellXGeneDownloadWorker(self, ids, name, self.parent.parent.moniker_recursively(name.replace(" ", "_") + "_corvo_RAW.h5ad", "datasets"))
        worker.signals.file_size.connect(disp_file_size)
        worker.signals.progress.connect(progress_bar)
        worker.signals.finished.connect(on_finished)
        # worker.signals.running.connect(on_running)
        self.cancel_download.connect(worker.stop)
        self.threadpool.start(worker)

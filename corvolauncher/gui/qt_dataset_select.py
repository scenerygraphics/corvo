import sys
import time
from functools import partial

from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal
from PyQt5.QtGui import QFont, QMouseEvent, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QWidget, QComboBox, QHBoxLayout, QPushButton, \
    QListView, QCompleter, QLineEdit, QFrame, QScrollArea

from corvolauncher.gui.job_runners.download_workers import DatasetDownloadWorker, DatasetInfoWorker
from corvolauncher.utilities.cellxgene_scrape import CellxGeneJSONRequests


class DatasetSelect(QWidget):
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        self.fetcher_layout = QVBoxLayout()
        self.fetcher_layout.setAlignment(Qt.AlignTop)
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
        # except Exception:
        except ValueError:  # experimental
            print("collection does not exist")

    @pyqtSlot(int)
    def load_collection_dsets(self, index):
        print("trying to load index: " + str(index))

        @pyqtSlot()
        def on_running():
            if self.fetcher_layout.count() > 3:  # on first loop
                self.fetcher_layout.removeWidget(self.children()[4])  # remove previous datasets
            self.fetcher_layout.addWidget(QLabel("fetching..."))

        @pyqtSlot(dict)
        def on_result(dset_map):
            print(dset_map)
            dset_button_layout = QVBoxLayout()
            dset_button_layout.setSpacing(15)
            dset_button_frame = QFrame()
            dset_button_layout.setAlignment(Qt.AlignTop)

            for dset in dset_map.keys():
                dset_layout = QHBoxLayout()
                dset_download_button = QPushButton("Download")

                dset_download_button.setObjectName(dset)
                # dset_download_button.setFixedWidth(65)
                # dset_download_button.setFixedHeight(20)
                dset_download_button.setMinimumSize(65, 25)
                dset_download_button.adjustSize()
                dset_download_button.clicked.connect(partial(self.launch_download, dset, dset_map))
                dset_layout.addWidget(dset_download_button)

                dset_cancel_button = QPushButton("Cancel")
                dset_cancel_button.setObjectName(dset + "cancel")
                dset_cancel_button.setMinimumSize(65, 25)
                dset_cancel_button.adjustSize()
                # dset_cancel_button.setFixedWidth(65)
                # dset_cancel_button.setFixedHeight(20)
                dset_cancel_button.hide()
                dset_layout.addWidget(dset_cancel_button)

                dset_label = QLabel(dset)
                dset_label.setMinimumHeight(25)
                dset_label.setFixedWidth(self.collection_container.width() - (dset_download_button.width() + 60))
                dset_label.adjustSize()
                dset_label.setWordWrap(True)
                dset_layout.addWidget(dset_label)
                dset_button_layout.addLayout(dset_layout)

            dset_button_layout.addStretch()
            # dset_button_layout.setMinimumSize(dset_button_layout.minimumSizeHint())
            dset_button_frame.setLayout(dset_button_layout)

            scroll_area = QScrollArea()
            scroll_area.setWidget(dset_button_frame)
            nb = len(dset_map.keys())
            if nb > 5:
                scroll_area.setFixedHeight(260)
            else:
                scroll_area.setFixedHeight((nb * 40) + 15)

            self.fetcher_layout.addWidget(scroll_area)

            self.fetcher_layout.removeWidget(self.children()[4])  # remove "fetching..." label

        if index != 0:  # in case hint is submitted as index
            worker = DatasetInfoWorker(self, index)
            worker.signals.running.connect(on_running)
            worker.signals.dict_result.connect(on_result)
            self.threadpool.start(worker)

    @pyqtSlot(str, dict)
    def launch_download(self, name, ids):
        print("object name:")
        print(self.findChild(QPushButton, name + "cancel"))

        @pyqtSlot()
        def on_finished():
            print("in finished")
            cancel_button.hide()
            download_button.show()
            cancel_button.setText("Cancel")  # reset button text
            self.parent.update_directory_index()

        @pyqtSlot(float)
        def disp_file_size(size):
            print(float(size) / 1000000)

        @pyqtSlot(int)
        def progress(done):
            sys.stdout.write("\r[%s%s]" % ('=' * done, ' ' * (50 - done)))
            sys.stdout.flush()

        @pyqtSlot()
        def cancel():
            # update the status bar from here
            worker.stop()
            cancel_button.disconnect()  # disconnect from temporary download worker as button is persistent
            cancel_button.setText("cancelling...")
            print("should be disconnected and only appear once")

        download_button = self.sender()
        cancel_button = self.findChild(QPushButton, name + "cancel")

        download_button.hide()
        cancel_button.show()
        cancel_button.clicked.connect(cancel)

        worker = DatasetDownloadWorker(self, name, ids)
        worker.signals.finished.connect(on_finished)
        worker.signals.file_size.connect(disp_file_size)
        worker.signals.progress.connect(progress)
        self.threadpool.start(worker)

    @pyqtSlot()
    def cancel(self):
        self.cancel_download.emit()

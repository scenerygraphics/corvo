import os
import sys
import time
import traceback
from functools import partial
from subprocess import Popen, DEVNULL, STDOUT

import requests
from PyQt5.QtCore import QRunnable, pyqtSlot, pyqtSignal, QObject, Qt
from PyQt5.QtWidgets import QVBoxLayout, QFrame, QHBoxLayout, QPushButton, QLabel

from corvolauncher.gui.job_runners.worker_signals import WorkerSignals
from corvolauncher.utilities.cellxgene_JSON_requests import CellxGeneJSONRequests


class CellXGeneDownloadWorker(QRunnable):

    def __init__(self, parent, ids, name, unique_name):
        super(CellXGeneDownloadWorker, self).__init__()

        self.parent = parent
        self.ids = ids
        self.name = name
        self.unique_name = unique_name

        self.shutdown = False
        self.signals = WorkerSignals()


        # only takes 0.3s to init and makes worker much simpler to instantiate each download

    @pyqtSlot()
    def run(self):
        try:
            self.signals.running.emit()
            # some datasets contain / in their name. Interpreted as directory and causes crash
            print("Downloading %s" % self.name)
            aws_hash = requests.post(
                "https://api.cellxgene.cziscience.com/dp/v1/datasets/" + self.ids[self.name]["dataset_id"] + "/asset/" +
                self.ids[self.name]["filetype_id"]).json()
            print(aws_hash)

            try:
                print(aws_hash["detail"])
                # can return json: {'detail': 'An internal server error has occurred. Please try again later.',
                # 'status': 500, 'title': 'Internal Server Error', 'type': 'about:blank'}
            except KeyError:
                h5ad_dataset = requests.get(aws_hash["presigned_url"], stream=True)
                # note some unusual characters cannot be parsed by the hdf5 reader

                # h5ad_dataset = requests.get(aws_hash["presigned_url"], stream=True)
                total_length = h5ad_dataset.headers.get('content-length')
                self.signals.file_size.emit(int(total_length))
                if total_length is None:  # no content length header
                    print("file has no contents")  # feed this back
                    pass
                    # f.write(h5ad_dataset.content)
                else:
                    with open("../resources/datasets/" + self.unique_name, "wb") as f:
                        dl = 0
                        total_length = int(total_length)
                        for data in h5ad_dataset.iter_content(chunk_size=round(total_length / 100)):
                            if not self.shutdown:
                                dl += len(data)
                                f.write(data)
                                done = int(50 * dl / total_length)
                                self.signals.progress.emit(done)
                            else:
                                os.remove("../resources/datasets/" + self.unique_name)
                                break

        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit()
        finally:
            self.signals.finished.emit()  # Done

    def stop(self):
        self.shutdown = True


class CellXGeneJSONWorker(QRunnable):
    def __init__(self, parent, index, *args, **kwargs):
        super(CellXGeneJSONWorker, self).__init__()

        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()
        self.parent = parent
        self.index = index

    @pyqtSlot()
    def run(self):
        if self.index != 0:
            dataset_map = self.parent.data_fetcher.get_collection_datasets(
                self.parent.collection_map[self.parent.collection_container.itemText(self.index)])

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
                dset_download_button.clicked.connect(partial(self.parent.launch_download, dset, dataset_map))
                dset_layout.addWidget(dset_download_button)

                # dset_cancel_button = QPushButton("Cancel")
                # dset_cancel_button.setObjectName(dset + "cancel")
                # dset_download_button.setFixedWidth(65)
                # dset_cancel_button.clicked.connect(self.cancel)
                # dset_cancel_button.hide()
                # dset_layout.addWidget(dset_cancel_button)

                dset_label = QLabel(dset)
                dset_label.setMinimumHeight(30)
                dset_label.setFixedWidth(self.parent.collection_container.width() - dset_download_button.width())
                # dset_label.setFixedWidth(self.collection_container.width() - 125)
                # dset_label.setMaximumWidth(self.collection_container.width() - dset_download_button.width())
                dset_label.setWordWrap(True)
                dset_layout.addWidget(dset_label)
                dset_button_layout.addLayout(dset_layout)
            dset_button_frame.setLayout(dset_button_layout)
            self.parent.fetcher_layout.addWidget(dset_button_frame)
        print("done")


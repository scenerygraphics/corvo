import os
import sys
import time
import traceback
from functools import partial
from subprocess import Popen, DEVNULL, STDOUT

import requests
from PyQt5.QtCore import QRunnable, pyqtSlot, pyqtSignal, QObject, Qt, QRect, QThread
from PyQt5.QtGui import QImage
from PyQt5.QtWidgets import QVBoxLayout, QFrame, QHBoxLayout, QPushButton, QLabel

from corvolauncher.gui.job_runners.worker_signals import WorkerSignals
from corvolauncher.utilities.cellxgene_scrape import CellxGeneJSONRequests


class DatasetDownloadWorker(QRunnable):
    """
    Attempts to download a dataset with a provided name and dictionary of corresponding IDs expected by the CellXGene
    AWS hash generator. Prints server error message if encountered (not yet linked to error slot).
        """

    def __init__(self, parent, name, ids):
        super(DatasetDownloadWorker, self).__init__()

        self.parent = parent
        self.ids = ids
        self.name = name

        self.unique_name = self.parent.parent.parent.moniker_recursively(self.name.replace(" ", "_") + "_corvo_RAW.h5ad", "datasets").replace("—", "-")

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
                # first checking that server error has not occurred
            except KeyError:
                # note some unusual characters cannot be parsed by the hdf5 reader
                h5ad_dataset = requests.get(aws_hash["presigned_url"], stream=True)
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
                                break
                    if self.shutdown:  # place here as file is used by open process in loop
                        os.remove("../resources/datasets/" + self.unique_name)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit()
        finally:
            print("emitting finished")
            self.signals.finished.emit()  # Done

    def stop(self):
        print("stopping")
        self.shutdown = True


class DatasetInfoWorker(QRunnable):
    """
    Fetch the dataset and filetype IDs for a given dataset located at a QComboBox index. Assumes the index and
    dataset names exist.
    """

    def __init__(self, parent, index, *args, **kwargs):
        super(DatasetInfoWorker, self).__init__()

        self.parent = parent
        self.index = index
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            self.signals.running.emit()
            dataset_map = self.parent.data_fetcher.get_collection_datasets(
                self.parent.collection_map[self.parent.collection_container.itemText(self.index)])
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.dict_result.emit(dataset_map)
        finally:
            self.signals.finished.emit()  # Done


class ModelDownloadWorker(QRunnable):
    def __init__(self):
        super(ModelDownloadWorker, self).__init__()

        pass
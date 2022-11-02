import sys
import traceback
import scanpy as sc

from PyQt5.QtCore import pyqtSlot, QObject, pyqtSignal, QRunnable
from corvolauncher.gui.job_runners.worker_signals import WorkerSignals


class AnndataWorker(QRunnable):
    """
    Constructs generic worker thread taking arguments and a function to run.
    """

    def __init__(self, dataset, *args, **kwargs):
        super(AnndataWorker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.dataset = dataset
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add callback to our kwargs
        # self.kwargs["progress_callback"] = self.signals.progress

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """
        # Retrieve args/kwargs here; and fire processing using them
        try:
            self.signals.running.emit()
            adata = sc.read_h5ad("../resources/datasets/" + self.dataset)
            obs = adata.obs
            var = adata.var

        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.adata_result.emit(obs, var)
        finally:
            self.signals.finished.emit()  # Done

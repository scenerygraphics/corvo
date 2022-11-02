import sys
import traceback

from PyQt5.QtCore import pyqtSlot, QObject, pyqtSignal, QRunnable

from corvolauncher.gui.job_runners.worker_signals import WorkerSignals
from corvolauncher.utilities.pre_process import PreProcess


class ProcessingWorker(QRunnable):
    """
    Constructs generic worker thread taking arguments and a function to run.
    """

    def __init__(self, dataset, *args, **kwargs):
        super(ProcessingWorker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.shutdown = False
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        self.dataset = dataset

        # Add the callback to our kwargs
        # self.kwargs["progress_callback"] = self.signals.progress

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """
        # Retrieve args/kwargs here; and fire processing using them
        try:
            self.signals.running.emit()
            PreProcess(self, self.dataset)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        finally:
            self.signals.finished.emit()  # Done

    def stop(self):
        self.shutdown = True

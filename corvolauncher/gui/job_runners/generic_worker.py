import sys
import traceback
from typing import Tuple

import requests
from PyQt5.QtCore import pyqtSlot, QObject, pyqtSignal, QRunnable
from corvolauncher.gui.job_runners.worker_signals import WorkerSignals


class GenericWorker(QRunnable):
    """
    Constructs generic worker thread taking arguments and a function to run.
    """

    def __init__(self, fn, *args, **kwargs):
        super(GenericWorker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
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
            # result = self.fn(*self.args, **self.kwargs)
            self.fn(*self.args, **self.kwargs)
        except Exception:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        # else:
        #     if type(result) == str:
        #         self.signals.str_result.emit(result)
        #     elif type(result) == requests.Response:
        #         self.signals.resp_result.emit(result)
        #     elif type(result) == int:
        #         self.signals.int_result.emit(result)
        #     else:
        #         pass
        finally:
            self.signals.finished.emit()  # Done

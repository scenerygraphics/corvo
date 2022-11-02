import pandas as pd
from PyQt5.QtCore import QObject, pyqtSignal, QMetaObject
from PyQt5.QtWidgets import QFrame, QVBoxLayout
from requests import Response


class WorkerSignals(QObject):
    """
    Worker signals
    """

    running = pyqtSignal()
    finished = pyqtSignal()
    cancelled = pyqtSignal()
    error = pyqtSignal()

    progress = pyqtSignal(int)
    file_size = pyqtSignal(int)

    result = pyqtSignal()
    adata_result = pyqtSignal(pd.DataFrame, pd.DataFrame)
    dict_result = pyqtSignal(dict)
    # int_result = pyqtSignal(int)
    # resp_result = pyqtSignal(Response)
    # str_result = pyqtSignal(str)

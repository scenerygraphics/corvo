from PyQt5.QtCore import QObject, pyqtSignal
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
    int_result = pyqtSignal(int)
    resp_result = pyqtSignal(Response)
    str_result = pyqtSignal(str)

from PyQt5.QtCore import QObject, pyqtSignal


class WorkerSignals(QObject):
    """
    Worker signals
    """

    running = pyqtSignal()
    finished = pyqtSignal()
    cancelled = pyqtSignal()
    error = pyqtSignal()
    # progress = pyqtSignal()

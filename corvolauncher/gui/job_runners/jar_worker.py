import time
from subprocess import Popen, DEVNULL, STDOUT

from PyQt5.QtCore import QRunnable, pyqtSlot, pyqtSignal, QObject

from corvolauncher.gui.job_runners.worker_signals import WorkerSignals


class JarWorker(QRunnable):
    def __init__(self, file_name, *args, **kwargs):
        super(JarWorker, self).__init__()

        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        self.file_name = file_name

        self.shutdown = False

    @pyqtSlot()
    def run(self):
        try:
            self.signals.running.emit()
            bash_command = "java -jar corvo-0.1.0-SNAPSHOT-all_3.jar " + "processed_datasets/" + self.file_name + \
                           " vosk-model-small-en-us-0.15"
            process = Popen(bash_command.split(), cwd="../resources")

            while not self.shutdown:
                time.sleep(0.1)

            process.kill()
            self.signals.cancelled.emit()

            while process.poll() is None:
                time.sleep(0.1)
        finally:
            self.signals.finished.emit()

    def stop(self):
        self.shutdown = True

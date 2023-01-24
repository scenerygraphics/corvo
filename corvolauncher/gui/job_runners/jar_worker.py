import os
import time
from pathlib import Path
from subprocess import Popen, DEVNULL, STDOUT
from importlib.resources import files

from PyQt5.QtCore import QRunnable, pyqtSlot, pyqtSignal, QObject

from corvolauncher.gui.job_runners.worker_signals import WorkerSignals


class JarWorker(QRunnable):
    def __init__(self, file_name, *args, **kwargs):
        super(JarWorker, self).__init__()

        self.jar_name = "corvo-0.1.0-SNAPSHOT-all.jar"
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        self.file_name = file_name

        self.shutdown = False

    @pyqtSlot()
    def run(self):
        try:
            self.signals.running.emit()

            base_path = os.path.join(str(Path.home()), ".corvo", "resources")
            jar_path = os.path.join(base_path, self.jar_name)
            if not os.path.exists(jar_path):
                print("Extracing scenery JAR, this is a one-time process...")
                jar_contents = files('corvolauncher').joinpath(self.jar_name).read_bytes()
                jar_file = open(jar_path, "wb")
                jar_file.write(jar_contents)
                jar_file.close()
                print("Extraction done, launching Corvo.")

            bash_command = "java -cp " + self.jar_name + " graphics.scenery.corvo.XVisualization processed_datasets/" + self.file_name + \
                           " vosk-model-small-en-us-0.15/vosk-model-small-en-us-0.15"
            print("Launching " + bash_command + " in directory " + base_path)
            process = Popen(bash_command.split(), cwd=base_path)

            while not self.shutdown:
                time.sleep(0.1)

            process.kill()
            # does not gracefully shut down, loading into jittery steamVR home screen
            self.signals.cancelled.emit()

            while process.poll() is None:
                time.sleep(0.1)
        finally:
            self.signals.finished.emit()

    def stop(self):
        self.shutdown = True

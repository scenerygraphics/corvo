import os.path
import subprocess

from corvolauncher.utilities import pre_process
from corvolauncher.utilities import jar_launch


class TaskManager:
    """
    main class to manage task order and pass data between class instances
    """

    def __init__(self, dataset: str):
        # 1: declare dataset paths for raw and processed
        self.dataset = dataset
        self.processed_dataset: str = "../resources/datasets/" + dataset.rstrip(".h5ad") + "_processed.h5ad"

        # 3: process dataset
        if not os.path.exists(self.processed_dataset):
            pre_process.PreProcess(self.dataset)

        # 4: launch Corvo with dataset
        bash_command = "java -jar corvo-0.1.0-SNAPSHOT-all.jar " + self.processed_dataset + \
                       " vosk-model-small-en-us-0.15"
        process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd="../resources")
        output, error = process.communicate()


if __name__ == "__main__":
    TaskManager("skin_MurisSenis.h5ad")

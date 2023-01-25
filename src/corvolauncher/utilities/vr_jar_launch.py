import os.path
import subprocess
from pathlib import Path


class JarLaunch:
    def __init__(self, data_name):
        self.path = os.path.join(str(Path.home()), ".corvo", "resources", "processed_datasets", data_name)
        print(self.path)
        bash_command = "java -jar corvo-0.1.0-SNAPSHOT-all.jar " + self.path + \
                       " vosk-model-small-en-us-0.15"
        process = subprocess.Popen(bash_command.split(), stdout=subprocess.PIPE, cwd="../resources")
        output, error = process.communicate()


if __name__ == "__main__":
    JarLaunch("marrow_processed.h5ad")

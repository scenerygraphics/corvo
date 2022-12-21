import json
import os.path
from pathlib import Path

import requests
from zipfile import ZipFile

from corvolauncher.utilities.speech_model_scrape import SpeechModelScrape


class Startup:
    """
    Check if Corvo has been run before by querying its system directories and if the speech model is present. If not
    found, directories are made and the model is downloaded and unzipped.
    """
    def __init__(self):
        # except NotADirectoryError:  # construct initial directories
        self.model_name = None
        self.speech_model_scrape = None

        if not os.path.isdir(os.path.join(str(Path.home()), ".corvo")):
            self.create_directories()

        if not os.path.isdir(os.path.join(str(Path.home()), ".corvo/resources/vosk-model-small-en-us-0.15")):
            self.download_speech_model()

    def create_directories(self):
        os.mkdir(os.path.join(str(Path.home()), ".corvo"))
        os.mkdir(os.path.join(str(Path.home()), ".corvo", "resources"))
        os.mkdir(os.path.join(str(Path.home()), ".corvo", "resources", "datasets"))
        os.mkdir(os.path.join(str(Path.home()), ".corvo", "resources", "processed_datasets"))

    def download_speech_model(self):
        self.speech_model_scrape = SpeechModelScrape()
        self.model_name = self.speech_model_scrape.data[0][1:][0][0]

        speech_model = requests.get("https://alphacephei.com/vosk/models/" + self.model_name + ".zip")
        with open(os.path.join(str(Path.home()), ".corvo", self.model_name + ".zip"), 'wb') as f:
            f.write(speech_model.content)

        with ZipFile(os.path.join(str(Path.home()), ".corvo", self.model_name + ".zip"), 'r') as zip_ref:
            zip_ref.extractall(os.path.join(str(Path.home()), ".corvo/resources", self.model_name))

        os.remove(os.path.join(str(Path.home()), ".corvo", self.model_name + ".zip"))


if __name__ == "__main__":
    Startup()

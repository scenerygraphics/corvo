from PyQt5.QtCore import pyqtSlot, Qt, pyqtSignal
from PyQt5.QtGui import QFont, QMouseEvent, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QVBoxLayout, QLabel, QWidget, QComboBox, QHBoxLayout, QPushButton, \
    QListView, QCompleter, QLineEdit, QFrame, QScrollArea

from corvolauncher.gui.job_runners.download_workers import DatasetDownloadWorker, DatasetInfoWorker
from corvolauncher.utilities.cellxgene_scrape import CellxGeneJSONRequests
from corvolauncher.utilities.speech_model_scrape import SpeechModelScrape


class SpeechModelSelect(QWidget):
    # implement check for if a model exists. If not, display the english model download buttons on launch and prompt
    # user to choose one. Until done, area should be tinted red.
    # Otherwise it is hidden and out of the way - can be expanded.
    # Change download layout to grid and include the model size.
    def __init__(self, parent, threadpool):
        super().__init__(parent)

        self.parent = parent
        self.threadpool = threadpool

        self.fetcher_layout = QVBoxLayout()
        self.fetcher_layout.setAlignment(Qt.AlignTop)
        self.setLayout(self.fetcher_layout)
        self.fetcher_layout.setSpacing(10)
        self.fetcher_layout.setAlignment(Qt.AlignLeft)

        #  give section titles a bold font
        self.title = QLabel("Vosk Speech Recognition Models")
        self.title_font = QFont()
        self.title_font.setPointSize(10)
        self.title_font.setBold(True)
        self.title.setFont(self.title_font)
        self.fetcher_layout.addWidget(self.title)

        self.speech_model_scrape = SpeechModelScrape()

        self.combo_font = QFont()
        self.combo_font.setPointSize(10)

        self.collection_container = QComboBox()
        self.collection_container.blockSignals(True)
        self.collection_container.currentIndexChanged.connect(self.load_models)
        self.collection_container.setFont(self.combo_font)
        self.collection_container.setFixedHeight(25)
        self.collection_container.setFixedWidth(self.parent.width())

        listView = QListView()
        listView.setWordWrap(True)
        listView.setSpacing(5)
        self.collection_container.setView(listView)
        self.collection_container.addItem("-- choose a model for speech recognition --")

        for language in self.speech_model_scrape.data:
            self.collection_container.addItem(language[0][0])
        self.fetcher_layout.addWidget(self.collection_container)

    @pyqtSlot(int)
    def load_models(self, index):
        print("trying to load index: " + str(index-1))
        if index != 0:
            l_list = self.speech_model_scrape.data[index - 1]
            language = l_list[0][0]
            language_models = l_list[1:]

            if self.fetcher_layout.count() > 2:  # on first loop
                self.fetcher_layout.removeWidget(self.children()[3])  # remove previous datasets

            model_button_layout = QVBoxLayout()
            model_button_layout.setSpacing(15)
            model_button_frame = QFrame()
            model_button_layout.setAlignment(Qt.AlignTop)

            for model in language_models:
                model_layout = QHBoxLayout()
                model_download_button = QPushButton("Download")

                model_download_button.setObjectName(model[0])
                model_download_button.setMinimumSize(65, 25)
                model_download_button.adjustSize()
                # model_download_button.clicked.connect(partial(self.launch_download, dset))
                model_layout.addWidget(model_download_button)

                cancel_download = QPushButton("Cancel")
                cancel_download.setObjectName(language + "cancel")
                cancel_download.setMinimumSize(65, 25)
                cancel_download.adjustSize()
                cancel_download.hide()
                model_layout.addWidget(cancel_download)

                model_label = QLabel(model[0])
                model_label.setMinimumHeight(25)
                model_label.setFixedWidth(self.collection_container.width() - (model_download_button.width() + 60))
                model_label.adjustSize()
                model_label.setWordWrap(True)
                model_layout.addWidget(model_label)
                model_button_layout.addLayout(model_layout)

            model_button_layout.addStretch()
            model_button_frame.setLayout(model_button_layout)
            self.fetcher_layout.addWidget(model_button_frame)
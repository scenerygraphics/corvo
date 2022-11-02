# CorvoLauncher
_CorvoLauncher_ is companion program for **Corvo** - a VR single-cell transcriptomics data explorer - that pre-processes ScanPy .h5ad datasets to be VR-ready and launches Corvo from a fat jar.

A Qt dataset manager and launcher is planned.

## Installation
It is recommended to use a virtual environment such as Anaconda for dependency management. This installation guide uses Anaconda to set up a virtual environment with the required dependencies. 

Install Anaconda here: https://www.anaconda.com/distribution/


### Clone Repository
```
git clone https://github.com/ljjh20/CorvoLauncher.git
```

### Create Conda Environment and Install Dependencies
```
conda create -n <name> python=3.7
conda activate <name>
pip install -r requirements.txt
```
### Vosk Language Model
You will need to download a Vosk English language model from: https://alphacephei.com/vosk/models and place it in **corvolauncher/resources**.

## Usage
### Launching
This prototype can be launched by running `corvolauncher/utilities/task_manager.py` with the dataset as an argument, and placing the raw dataset into `/corvolauncher/resources/datasets`.

### Finding Datasets
A large repository of compatible single-cell RNA transcriptomics datasets are available from the Chan Zuckerberg Initiative (must be downloaded in .h5ad format): 
https://cellxgene.cziscience.com/datasets


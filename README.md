# Corvo
_Corvo_ is a free and open-source tool that takes the visualization and analysis of CZI CellxGene (https://cellxgene.cziscience.com/) single-cell datasets to 3D Virtual Reality. Similar to CellxGene, Corvo takes a no-code approach for the end user, but also offers multimodal user input to facilitate fast navigation and analysis, and is interoperable with the existing Python data science ecosystem

You can launch Corvo by either running ```corvolauncher/gui/qt_main_window.py```, or by launching the accompanying app that is installed when you run ```setup.py```. 

## Installation
Using an environment manager is highly reccommended. The following instructions use Anaconda, but you are welcome to use another. Pip is also required for installing dependencies.
### Clone Repository
```
https://github.com/scenerygraphics/corvo.git
```
### Run setup.py found in the root directory
```
python setup.py
```
### Create Conda Environment and Install Dependencies
```
conda create -n <name> python=3.7
conda activate <name>
pip install -r requirements.txt
```

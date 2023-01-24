# Corvo
_Corvo_ is a free and open-source tool that takes the visualization and analysis of [CZI CellxGene](https://cellxgene.cziscience.com/) single-cell datasets to 3D Virtual Reality. Similar to CellxGene, Corvo takes a no-code approach for the end user, but also offers multimodal user input to facilitate fast navigation and analysis, and is interoperable with the existing Python data science ecosystem

![Corvo’s VR interface: A: menu (1), press-
to-talk voice interface, with genes recognized from voice input
listed (2), list of precomputed genes (3) B: UMAP plot with
gene expression scale (1) C: interaction / selection sphere D:
annotation key.](./images/corvo.png)
Examples of Corvo’s VR interface: A: menu (1), press-
to-talk voice interface, with genes recognized from voice input
listed (2), list of precomputed genes (3) B: UMAP plot with
gene expression scale (1) C: interaction / selection sphere D:
annotation key.


You can launch Corvo by either running ```corvolauncher/gui/qt_main_window.py```, or by launching the accompanying app that is installed when you run ```setup.py```. 

Corvo is currently on [arXiv](https://arxiv.org/abs/2212.00519). Please cite this preprint in case you are using Corvo.

## Requirements

* A SteamVR-supported VR headset (e.g. Oculus Quest 2, HTC Vive) - Corvo not work without a VR headset at the moment.
* A computer running Windows 10, equipped with a Geforce GTX 1050-class GPU at minimum, and 8 GB RAM. Corvo might
  also work with SteamVR on Linux, this is however untested at the moment.

## Installation from source code
The following instructions use [Anaconda](https://www.anaconda.com/), but you are welcome to use an alternative, like miniconda or mamba.

* Clone repository, recursive cloning here is important, as this will also clone the [corvo-core](https://github.com/scenerygraphics/corvo-core) submodule.
  ```shell
  git clone --recursive https://github.com/scenerygraphics/corvo
  ```
* Import the provided `environment.yml` file for conda. This will install the required dependencies and create a new environment named `corvo`.
  ```shell
  conda env create -f environment.yml
  conda activate corvo
  ```
* Run the provided setup script. This will compile `corvo-core`, and use setuptools to install corvo.
  ```shell
  # On Linux and macOS
  ./build.sh
  # on Windows
  .\bld.bat
  ```
* Finally, you can run Corvo, by issueing `corvo` on the command line. Enjoy 👍

## Installation via conda (experimental)

Corvo can also be installed via conda directly. This skips the need for the build process outlined above. Run the following commands:
```shell
conda install -c skalarproduktraum corvo
```

Afterwards, you can run Corvo by issueing `corvo` on the command line.

## Conda package creation (experimental, optional)

You can also create your own custom conda package, by running the following command in the directory of this repository:
```shell
conda-build .
```
This requires `conda-build` to be installed and will build a conda package for the architecture you're currently running on.

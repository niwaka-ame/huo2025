# Code of figures in Huo _et al._ _Commun. Biol._ (2025)
This repository is a step-by-step guide to reproduce the figures in Huo _et al._ _Commun. Biol._ (2025).

## Prerequisites
1. Operating system: Linux, Mac OS or [Windows Subsystem for Linux 2 (WSL2)](https://learn.microsoft.com/en-us/windows/wsl/install).
2. Anaconda or Miniconda. For example, please refer to [the official documentation of miniconda](https://docs.anaconda.com/miniconda/install/).
3. `pipx` and `poetry`. Please refer to [the official installation guide of Poetry](https://python-poetry.org/docs/#installation)

As an example, below are the commands for installing `pipx` and `poetry` on Fedora Linux 41:
```sh
sudo dnf install pipx
pipx ensurepath
pipx install poetry
```

## Installing packages
1. Type the following commands one by one in your Terminal:
```sh
conda create -n huo2025fig python=3.10
conda activate huo2025fig
cd ~/Downloads/huo2025  # change this to where the codes downloaded from Zenodo are stored
git submodule update --init
poetry install
```
This environment will suffice to reproduce all figures except S2, S6, S7 and S11.

2. To reproduce Fig S2, please temporarily downgrade the `pandas` package (required for an old version of the `wela` package):
```sh
pip install -U pandas==1.4.4
```
To upgrade `pandas` back to the version specified by `poetry` after reproducing Fig S2:
```sh
cd ~/Downloads/huo2025  # make sure we're in the codes' folder
poetry install
```

3. To reproduce Fig S6, S7 and S11, we need to install R, RStudio and the following packages: `tidyverse`, `rstudioapi`, `BiocManager`, `devtools`, `repr`, `ggplot2`, `DESeq2`, `tidyqpcr`. 

## Reproduce the figures
1. Now under the folder `~/huo2025/` there should be three folders: `data`, `src` and `svg`.
2. Create a new folder to store the reproduced figures:
``` sh
cd ~/Downloads/huo2025 # where the codes are stored
mkdir fig
```
3. Run the Python files one by one in the `src` folder, e.g.
```sh
cd ~/Downloads/huo2025 # where the codes are stored
python src/fig1.py
```
4. To reproduce Figures S6 and S7, make sure you run the R codes first in the `src/intermediate` folder. The output panels will be stored in the `fig` folder.

## Note
The code for Figure S10 will send request to the YeastEnrichR server; please pay attention to the frequency of request. 

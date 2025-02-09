# Code for figures in Huo _et al._ (2025). The type of carbon source not the growth rate it supports can determine diauxie in _Saccharomyces cerevisiae_. _Communications Biology_. 
This repository is a step-by-step guide to reproduce the figures in Huo _et al._ _Commun. Biol._ (2025).

## Prerequisites
1. Operating system: Linux, Mac OS or [Windows Subsystem for Linux 2 (WSL2)](https://learn.microsoft.com/en-us/windows/wsl/install).
2. [Git](https://git-scm.com/). This usually ships with the above systems.
3. Anaconda or Miniconda. For example, please refer to [the official documentation of miniconda](https://docs.anaconda.com/miniconda/install/).
4. `pipx` and `poetry`. Please refer to [the official installation guide of Poetry](https://python-poetry.org/docs/#installation)

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
mkdir ~/temp  # where you want to temporarily store this code repository
cd ~/temp
git clone https://git.ecdf.ed.ac.uk/s1856140/huo2025.git
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
cd ~/temp/huo2025  # make sure we are in the code repository
poetry install
```

3. To reproduce Fig S6, S7 and S11, we need to install R, RStudio and the following packages: `tidyverse`, `rstudioapi`, `BiocManager`, `devtools`, `repr`, `ggplot2`, `DESeq2`, `tidyqpcr`. 

## Fetching the data from Zenodo
The data needed to reproduce the figures are available on Zenodo: https://doi.org/10.5281/zenodo.14840343. Download and move `data.zip` into the code repository, e.g. `~/temp/huo2025`, and then decompress it. 

For example, if you use the Terminal to do so:
``` sh
mv ~/Downloads/data.zip ~/temp/huo2025 
cd ~/temp/huo2025
unzip data.zip  # you may need to install `unzip` if missing 
```

## Reproducing the figures
1. Now under the folder there should be four folders: `data`, `src`, `fig` and `svg`.

2. Run the Python files one by one in the `src` folder, e.g.
```sh
cd ~/temp/huo2025 # the code repository
python src/fig1.py
```
The figure will be output to the `fig` folder.

3. To reproduce Figures S6 and S7, make sure you run the R code first in the `src/intermediate` folder. The panels will be output to the `fig` folder.

## Note
The code for Figure S10 will send request to the [YeastEnrichR server](https://maayanlab.cloud/YeastEnrichr/); please pay attention to the frequency of request. 

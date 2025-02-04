# Code of figures in Huo _et al._ _Commun. Biol._ (2025)
This repository is a step-by-step guide to reproduce the figures in Huo _et al._ _Commun. Biol._ (2025).

## Prerequisites
1. Operating system: Linux, Mac OS or Windows Subsystem for Linux 2 (WSL2).
2. Anaconda or Miniconda

## Installing packages
1. Type the following commands one by one in your Terminal:
```sh
cd ~
git clone https://git.ecdf.ed.ac.uk/s1856140/huo2025.git
conda create -n huo2025 python=3.10
conda activate huo2025
pip install poetry
cd huo2025
git submodule update --init
poetry install
```
This environment will suffice to reproduce all figures except S2, S6, S7 and S11.

2. To reproduce Fig S2, please temporarily downgrade the `pandas` package (which is required for an old version of the `wela` package):
```sh
pip install -U pandas==1.4.4
```

3. To reproduce Fig S6, S7 and S11, we need to install R with the following packages: `tidyverse`, `ggplot2`, `BiocManager`, `DESeq2`, `devtools`, `tidyqpcr`, `repr`. 

## Fetch the data 
1. Download Supplementary Data 2.
2. Unzip and move to `~/huo2025/`.
3. Rename the folder as `data`.

## Reproduce the figures
1. Now under the folder `~/huo2025/` there should be three folders: `data`, `src` and `svg`.
2. Create a new folder to store the reproduced figures: `mkdir ~/huo2025/fig`
3. Run the Python files one by one in the `src` folder, e.g.
```sh
cd ~/huo2025/src
python fig1.py
```
4. To reproduce Figures S6 and S7, make sure you run the R codes first in the `src/intermediate` folder.

## Note
The code for Figure S10 will send request to the YeastEnrichR server; please pay attention to the frequency of request. 

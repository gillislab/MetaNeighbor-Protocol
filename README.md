# Scaling up reproducible research for single cell transcriptomics using MetaNeighbor (Protocols)

In this repository, we present 3 protocols to illustrate how MetaNeighbor can be used to quantify and characterize cell type replicability across single cell transcriptomic datasets.

## Installation and data download
To run the protocols, you will need to install the R or Python version of MetaNeighbor. Installation instructions are included in Protocol 1, but you can also install the development version of MetaNeighbor for [R](https://github.com/gillislab/MetaNeighbor) or [Python](https://github.com/gillislab/pyMN) by following the corresponding links (installation instructions are included in the README). All protocols have been tested on Windows 10, MacOS Catalina 10.15 and Linux RHEL7 and are expected to run on an up-to-date R (tested on versions 3.6 and 4.0) or Python environment.

This repository also includes the scripts generating the data used in the protocols (data/ directory). Alternatively, you can access ready-to-use data from FigShare for [R](https://doi.org/10.6084/m9.figshare.13020569) or [Python](https://doi.org/10.6084/m9.figshare.13034171).

## Protocol description

We provide 3 protocols and a procedure for troubleshooting common issues ("anticipated results") in various notebook formats:
 - In Protocol 1, we show how to use MetaNeighibor to quantify the replicability of cell types across datasets. MetaNeighbor quantifies the robustness of cell types with respect to various sources of technical variability, such as lab of origin, sequencing technology or clustering pipeline. The quantification enables to differentiate between robust cell types (easily detected, independent of experimental conditions) and cell types that require more careful characterization. 
 - In Protocol 2, we show how to use MetaNeighbor to rapidly assess newly annotated cell types against a large reference taxonomy. Protocol 2 has similar aims as Protocol 1, but considerably speeds up computations by pre-training a model on the reference taxonomy. Specifically, we show how neuronal cell types can be assessed within seconds against a reference taxonomy from the Brain Initiative Cell Census Network (BICCN) containing hundreds of cell types determined from half a million high quality cells.
 - In Protocol 3, we show how to use MetaNeighbor to functionally characterize cell type replicability. Once a replicating cell type has been identified (in Protocol 1 or 2), we can use MetaNeighbor’s supervised framework to identify which pathways contribute most to the cell type’s identity. Within one hour, MetaNeighbor can process thousands of functional gene sets and highlight the most relevant pathways.
 - In Anticipated Results, we show how to diagnose and troubleshoot common issues, such as problems related to a bad selection of highly variable genes (the only "parameter" in MetaNeighbor).
 
For each protocol, we propose the following formats:
 - Rmd: Rmarkdown notebook format, compatible with Rstudio and Jupyter with R kernel.
 - ipynb (R): Jupyter notebook format (for Jupyter with R kernel).
 - ipynb (python, located in python_notebooks directoly): Jupyter notebook format (for Jupyter with python kernel).
 - PDF: rendered version of notebooks, including output of the different code blocks.
 
The expected run time is as follows:
 - Installation (Protocol 1, Step 0): usually 1-2 minutes, but can go up to 20 minutes when starting from a completely empty R environment.
 - Protocol 1: 1-5 minutes.
 - Protocol 2: 1-5 minutes.
 - Protocol 3: 30-90 minutes.

Note that Protocol 3 is a little more memory intensive than Protocols 1 and 2. 8GB are sufficient to run the first two protocols, but we have found that 16GB may be necessary to run Protocol 3.

# Scaling up reproducible research for single cell transcriptomics using MetaNeighbor (Protocol)

In this repository, we present 3 procedures to illustrate how MetaNeighbor can be used to quantify and characterize cell type replicability across single cell transcriptomic datasets.

## Installation and data download
To run the protocol, you will need to install the R or Python version of MetaNeighbor. Installation instructions are included in the [Equipment Setup](equipment_setup.Rmd), but you can also install the development version of MetaNeighbor for [R](https://github.com/gillislab/MetaNeighbor) or [Python](https://github.com/gillislab/pyMN) by following the corresponding links (installation instructions are included in the README). All proceduress have been tested on Windows 10, MacOS Catalina 10.15 and Linux RHEL7 and are expected to run on a reasonably up-to-date R (tested on versions 3.6 and 4.0) or Python environment.

The [data/](data/) repository includes the scripts generating the data used in the protocol. Alternatively, you can access ready-to-use data from FigShare for [R](https://doi.org/10.6084/m9.figshare.13020569) or [Python](https://doi.org/10.6084/m9.figshare.13034171).

The [pretrained_models/](pretrained_models/) repository includes ready-to-use models of the BICCN taxonomy for Procedure 2 ("Assessing cell type replicability against a pre-trained reference taxonomy"). Procedure 2 explains how to create these models in Section 1 ("Pre-train a reference MetaNeighbor model"), but you can also download the pre-trained models and skip to Section 2 ("Compare annotations to pre-trained taxonomy").

The [anticipated_results_figs/](anticipated_results_figs/) repository includes scripts used to estimate the effect of batch effects on MetaNeighbor results in three scenarios: lower sensitivity (lower number of UMIs per cell), higher noise, and multi-modal analyses.


## Protocol description

We provide an installation procedure, 3 main procedures and a procedure for troubleshooting common issues ("anticipated results") in various notebook formats:
 - In Equipment Setup, we show how to install MetaNeighbor and other packages necessary to run the protocol.
 - In Procedure 1, we show how to use MetaNeighbor to quantify the replicability of cell types across datasets. MetaNeighbor quantifies the robustness of cell types with respect to various sources of technical variability, such as lab of origin, sequencing technology or clustering pipeline. The quantification enables to differentiate between robust cell types (easily detected, independent of experimental conditions) and cell types that require more careful characterization. 
 - In Procedure 2, we show how to use MetaNeighbor to rapidly assess newly annotated cell types against a large reference taxonomy. Procedure 2 has similar aims as Procedure 1, but considerably speeds up computations by pre-training a model on the reference taxonomy. Specifically, we show how neuronal cell types can be assessed within seconds against a reference taxonomy from the Brain Initiative Cell Census Network (BICCN) containing hundreds of cell types determined from half a million high quality cells.
 - In Procedure 3, we show how to use MetaNeighbor to functionally characterize cell type replicability. Once a replicating cell type has been identified (in Procedure 1 or 2), we can use MetaNeighbor’s supervised framework to identify which pathways contribute most to the cell type’s identity. Within one hour, MetaNeighbor can process thousands of functional gene sets and highlight the most relevant pathways.
 - In Anticipated Results, we show how to diagnose and troubleshoot common issues, such as problems related to a bad selection of highly variable genes (the only "parameter" in MetaNeighbor).
 
For each procedure, we propose the following formats:
 - Rmd: Rmarkdown notebook format, compatible with Rstudio and Jupyter with R kernel.
 - ipynb (R): Jupyter notebook format (for Jupyter with R kernel).
 - ipynb (python, located in python_notebooks directoly): Jupyter notebook format (for Jupyter with python kernel).
 - PDF: rendered version of notebooks, including output of the different code blocks.
 
The expected run time is as follows:
 - Equipment Setup (Installation): usually 1-2 minutes, but can go up to 20 minutes when starting from a completely empty R environment.
 - Procedure 1: 1-5 minutes.
 - Procedure 2: 1-5 minutes.
 - Procedure 3: 30-90 minutes.

Note that Procedure 3 is a little more memory intensive than Procedures 1 and 2. 8GB are sufficient to run the first two procedures, but we have found that 16GB may be necessary to run Procedure 3.

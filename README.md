# Scaling up reproducible research for single cell transcriptomics using MetaNeighbor (Protocols)

In this repository, we present 3 protocols to illustrate how MetaNeighbor can be used to quantify and characterize cell type replicability across single cell transcriptomic datasets.

To run the protocols, you will need to install the R or Python version of MetaNeighbor. Installation instructions are included in Protocol 1, but you can also install the development version of MetaNeighbor for [R](https://github.com/gillislab/MetaNeighbor) or [Python](https://github.com/gillislab/pyMN) by following the corresponding links (installation instructions are included in the README). All protocols have been tested on Windows 10, MacOS Catalina 10.15 and Linux RHEL7 and are expected to run on a reasonably up-to-date R or Python environment.

The [data/] repository includes the scripts generating the data used in the protocols. Alternatively, you can access ready-to-use data from FigShare for [R](https://doi.org/10.6084/m9.figshare.13020569) or [Python](https://doi.org/10.6084/m9.figshare.13034171).

The [pretrained_models/] repository includes ready-to-use models of the BICCN taxonomy for Protocol 2 ("Assessing cell type replicability against a pre-trained reference taxonomy"). Protocol 2 explains how to create these models in Step 1 ("Pre-train a reference MetaNeighbor model"), but you can also download the pre-trained models and skip to Step 2 ("Compare annotations to pre-trained taxonomy").


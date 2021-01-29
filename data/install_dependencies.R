
if (!require('devtools')) {
  install.packages('devtools', quiet=TRUE)
}
devtools::install_github("gillislab/MetaNeighbor")

to_install = c("SingleCellExperiment", "Matrix", "rhdf5", "tidyverse",
               "org.Mm.eg.db", "GO.db")
installed = sapply(to_install, requireNamespace, quietly = TRUE)
if (sum(!installed) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", quiet = TRUE)
        BiocManager::install()
    }
    BiocManager::install(to_install[!installed])
}

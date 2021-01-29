library(rhdf5)
library(Matrix)
library(tidyverse)
library(SingleCellExperiment)


main = function() {
    saveRDS(make_atac_dataset(), "snatac/snatac.rds")
    saveRDS(make_snmc_dataset(), "snmc/snmc.rds")
    saveRDS(make_multi_modal_dataset(), "biccn_multi_modal_hvg.rds")
}

make_atac_dataset = function(atac_dir = "snatac") {
    metadata = read_csv(file.path(atac_dir, "metadata.csv.gz")) %>%
        mutate(cell_name = paste(sample, barcode, sep = "."))
    sample_name = unique(metadata$sample)
    result = lapply(sample_name, function(s) {
        read_counts_from_snap(file.path(atac_dir, paste0(s, ".snap")), s)
    })    
    result = do.call(cbind, result)
    result = result[, metadata$cell_name]
    
    result = SingleCellExperiment(list(counts = result))
    result$cluster_label = metadata$SubCluster
    result$subclass_label = metadata$MajorCluster
    result$class_label = atac_subclass_to_class(result$subclass_label)
    
    return(result)
}

read_counts_from_snap = function(filename, sample_name) {
    counts = h5read(filename, "GM")
    cell_names = paste(sample_name, h5read(filename, "BD/name"), sep = ".")
    result = Matrix::sparseMatrix(i = as.integer(counts$idy),
                                  j = as.integer(counts$idx),
                                  x = as.integer(counts$count),
                                  dimnames = list(as.character(counts$name), cell_names))
    return(result)
}

atac_subclass_to_class = function(original_name) {
    map = c()
    map[c("OGC", "Smc", "OPC", "MGC", "Endo", "ASC")] = "Non-Neuronal"
    map[c("CGE", "Sst", "Pv")] = "GABAergic"
    map[c("L6.CT", "L4", "L23.c", "L5.IT.b", "NP", "L23.a", "L5.IT.a", "L5.PT", "L6.IT", "L23.b")] = "Glutamatergic"
    map["Other"] = "Other"
    result = map[original_name]
}

make_snmc_dataset = function(snmc_dir = "snmc") {
    metadata = read_csv(file.path(snmc_dir, "metadata.csv.gz"))
    sample_name = unique(paste(metadata$Region, metadata$FACS_Date, sep="-"))
    result = lapply(sample_name, function(s) {
        read_counts_from_mcds(file.path(snmc_dir, paste0(s, ".mcds")))
    })    
    result = do.call(cbind, result)    
    result = result[, metadata$index]
    rownames(result) = as.character(rownames(result))
    
    result = SingleCellExperiment(list(counts = result))
    result$cluster_label = metadata$SubCluster
    result$subclass_label = metadata$MajorCluster
    result$class_label = snmc_subclass_to_class(result$subclass_label)
    return(result)
}

read_counts_from_mcds = function(filename) {
    counts = h5read(filename, "/")
    # dimension 1: mc_type 1="mc", 2="cov" 
    # dimension 2: strand_type 1="both"
    # dimension 3: count_type 1="CGN", 2="CHN"
    mc = counts$gene_da[1,1,2,,]
    cov = counts$gene_da[2,1,2,,]
    result = Matrix::Matrix(mc / (cov+1), sparse = TRUE)
    dimnames(result) = list(counts$gene, counts$cell)
    return(result)
}

snmc_subclass_to_class = function(original_name) {
    result = original_name
    result[startsWith(original_name, "L")] = "Glutamatergic"
    result[startsWith(original_name, "CGE")] = "GABAergic"
    result[startsWith(original_name, "MGE")] = "GABAergic"
    result[startsWith(original_name, "NonN")] = "Non-Neuronal"
    return(result)
}

make_multi_modal_dataset = function(){
    rna = readRDS("full_biccn_hvg.rds")

    snmc = readRDS("snmc/snmc.rds")
    counts(snmc) = -counts(snmc)
    snmc$study_id = "snmc"
    gene_symbols = convert_gene_names(rownames(snmc))
    is_na = is.na(gene_symbols)
    snmc = snmc[!is_na,]
    rownames(snmc) = gene_symbols[!is_na]
    
    atac = readRDS("snatac/snatac.rds")
    atac$study_id = "atac"
    
    full_data = list(rna, atac, snmc)
    full_data = MetaNeighbor::mergeSCE(full_data)
    full_data$study_id = c(rna$study_id, atac$study_id, snmc$study_id)
    return(full_data)
}   

convert_gene_names = function(snmc_names) {
    gene_annotation = read.table("gene_annotation/mgi_ensembl_191216.txt",
                                 sep = "\t", quote = "", comment.char = "",
                                 stringsAsFactors = FALSE, header = FALSE)
    gene_map = gene_annotation[,2]
    names(gene_map) = gene_annotation[,6]
    snmc_ensembl = sapply(strsplit(snmc_names, split = ".", fixed = TRUE), head, 1)
    result = gene_map[snmc_ensembl]
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
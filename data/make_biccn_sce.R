
library(SingleCellExperiment)
library(Matrix)
library(rhdf5)
library(tidyverse)
library(MetaNeighbor)


main = function() {
    #make_all_datasets()
    full_dataset = make_merged_dataset()
    hvgs = make_hvgs(full_dataset)
    make_hvg_dataset(full_dataset, hvgs)
    make_gaba_dataset(full_dataset)
}

make_all_datasets = function(output_dir = ".") {
    dataset_names = c("scCv2", "scCv3", "scSS", "snCv2", "snCv3M", "snCv3Z", "snSS")
    dir.create(output_dir, showWarnings=FALSE)
    for (input_dir in dataset_names) {
        print(paste0("Creating dataset ", input_dir, "..."))
        counts = if (is_SS(input_dir)) {
            load_smart_counts(input_dir)
        } else {
            load_10x_counts(input_dir)
        }
        metadata = load_metadata(input_dir) %>%
            inner_join(load_clusters(input_dir), by = "sample_id", suffix = c("_from_metadata", "")) %>%
            filter(passed_qc == TRUE) %>%
            left_join(load_joint_clusters(input_dir), by = "sample_id")
        # in some of the 10X datasets, sample names in data and metadata do not match exactly
        # in these cases, we need to remove the library name from the sample name
        sample_ids = metadata$sample_id
        if (!any(sample_ids %in% colnames(counts))) {
            sample_ids = metadata %>%
                rowwise() %>%
                transmute(short_id = gsub(Lib_Name, "", sample_id, fixed = TRUE)) %>%
                pull()
        }
        counts = counts[, sample_ids]
        
        sce = SingleCellExperiment(list(counts = counts), colData = metadata)
        sce$study_id = input_dir
        saveRDS(sce, file.path(input_dir, paste0(input_dir, ".rds")))
        rm(sce); rm(counts); rm(metadata); gc();
    }
}

is_SS = function(dataset_name) {
    grepl("SS", dataset_name)
}

load_10x_counts = function(subdir) {
  h5_data = h5read(file.path(subdir, "umi_counts.h5"), "/")[[1]]
  gene_names = as.character(h5_data$features$name)
  barcodes = as.character(h5_data$barcodes)
  result = Matrix::sparseMatrix(i = h5_data$indices+1,
                                p = h5_data$indptr,
                                x = h5_data$data,
                                dim = c(length(gene_names), length(barcodes)),
                                dimnames = list(gene_names, barcodes))
  return(result)
}

load_smart_counts = function(subdir) {
  result = read.table(file.path(subdir, "exon.counts.csv.gz"),
                      header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  return(Matrix(as.matrix(result), sparse = TRUE))
}

load_metadata = function(subdir) {
    filename = if (is_SS(subdir)) {"sample_metadata.csv.gz"} else {"sample_metadata.csv"}
    result = read.csv(file.path(subdir, filename), stringsAsFactors = FALSE, row.names=1)
    if (subdir != "snCv3M") {
        result = rownames_to_column(result, "sample_id")
    } else {
        result$sample_id = result$sample_name
    }
    return(result)
}

load_clusters = function(subdir, check_qc=TRUE) {
    clusters = read.csv(file.path(subdir, "cluster.membership.csv"), stringsAsFactors = FALSE)
    colnames(clusters) = c("sample_id", "cluster_id")
    cluster_info = read.csv(file.path(subdir, "cluster.annotation.csv"), stringsAsFactors = FALSE)
    result = inner_join(clusters, cluster_info, by="cluster_id")
    if (check_qc) {
        result$passed_qc = result$sample_id %in% qc_cells(subdir)
    }
    return(result)
}

qc_cells = function(subdir) {
    qc_cells = read.csv(file.path(subdir, "QC.csv"), stringsAsFactors = FALSE)$x
    return(qc_cells)
}

load_joint_clusters = function(dataset_name) {
    result = load_clusters("joint_annotation", check_qc = FALSE) %>%
        dplyr::rename_all(function(x) {paste0("joint_", x)}) %>%
        separate(joint_sample_id, c("dataset", "sample_id"), sep = "\\.") %>%
        mutate(dataset = convert_dataset_name(dataset)) %>%
        filter(dataset == dataset_name) %>%
        select(-dataset)
    return(result)
}

convert_dataset_name = function(dataset_name) {
    new_names = c("10X_cells_v2_AIBS"="scCv2", "10X_cells_v3_AIBS"="scCv3", "10X_nuclei_v2_AIBS"="snCv2",
                  "10X_nuclei_v3_AIBS"="snCv3Z", "10X_nuclei_v3_Broad"="snCv3M",
                  "SmartSeq_cells_AIBS"="scSS", "SmartSeq_nuclei_AIBS"="snSS")
    return(new_names[dataset_name])
}

make_merged_dataset = function() {
    print("Creating full dataset...")
    dataset_names = c("scCv2", "scCv3", "scSS", "snCv2", "snCv3M", "snCv3Z", "snSS")
    dataset = lapply(set_names(dataset_names), function(n) { readRDS(file.path(n, paste0(n, ".rds"))) })
    gc()
    dataset = mergeSCE(dataset)
    gc()
    saveRDS(dataset, "full_biccn.rds")
    return(dataset)
}

make_hvgs = function(full_dataset) {
    # to always obtain the same set of HVGs despite downsampling,
    # we choose an arbitrary random seed
    set.seed(17)
    system.time({hvg = variableGenes(full_dataset, exp_labels = full_dataset$study_id, downsampling_size = 10000)})
    write(hvg, "biccn_hvgs.txt")
    return(hvg)
}

make_hvg_dataset = function(full_dataset, hvg) {
    print("Creating HVG-restricted dataset...")
    saveRDS(full_dataset[hvg,], "full_biccn_hvg.rds")
}

make_gaba_dataset = function(full_dataset) {
    print("Creating GABAergic neuron-restricted dataset...")
    dataset = full_dataset[, full_dataset$joint_class_label == "GABAergic"]
    saveRDS(dataset, "biccn_gaba.rds")
}

main()
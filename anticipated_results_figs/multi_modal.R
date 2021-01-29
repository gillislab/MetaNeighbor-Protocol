
library(tidyverse)
library(SingleCellExperiment)
library(MetaNeighbor)


main = function() {
    pdf("multi_modal.pdf")
    make_multi_modal_figures()
    dev.off()
}

make_multi_modal_figures = function() {
    full_data = readRDS("../data/biccn_multi_modal_hvg.rds")
    full_data = full_data[, full_data$class_label %in% c("GABAergic", "Glutamatergic", "Non-Neuronal")]

    # Basic MN analysis
    hvg = rownames(full_data)
    auroc = MetaNeighborUS(hvg, full_data, study_id = full_data$study_id,
                           cell_type = full_data$subclass_label,
                           fast_version = TRUE)

    base_cols = RColorBrewer::brewer.pal(n = 3, "Set1")
    names(base_cols) = c("atac", "methylation", "rna")
    dataset_to_col = rep(base_cols["rna"], length(unique(full_data$study_id)))
    names(dataset_to_col) = unique(full_data$study_id)
    dataset_to_col["atac"] = base_cols["atac"]
    dataset_to_col["snmc"] = base_cols["methylation"]
    plotHeatmap(auroc,
                ColSideColors = dataset_to_col[getStudyId(colnames(auroc))],
                RowSideColors = dataset_to_col[getStudyId(colnames(auroc))])
    
    # remove same-study similarities
    #study_matrix = matrix(getStudyId(rownames(auroc)), nrow = nrow(auroc), ncol = ncol(auroc))
    #is_same_study = study_matrix == t(study_matrix)
    #auroc_na = auroc
    #auroc_na[is_same_study] = NA
    #plotHeatmap(auroc_na)
    #ggPlotHeatmap(auroc_na, label_size = 5)
    
    # plot # best hits / AUROC for each study
    best_hits = find_best_hits(auroc) %>%
        mutate(ref_study = getStudyId(ref_cell_type))
    
    to_plot = best_hits %>%
        group_by(ref_study, ref_cell_type) %>%
        summarize(n_hits = sum(is_reciprocal)) %>%
        group_by(ref_study) %>% 
        summarize(m = mean(n_hits), sd = sd(n_hits)) %>%
        mutate(data_type = data_type(ref_study)) %>%
        ggplot(aes(x = ref_study, y = m, fill = data_type)) +
        geom_col() +
        geom_linerange(aes(ymin = m-sd, ymax = m + sd)) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        labs(x = "Reference dataset", y = "Number of reciprocal best hits",
             fill = "Data type") +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)

    to_plot = best_hits %>%
        group_by(ref_study, ref_cell_type) %>%
        summarize(auroc = mean(auroc)) %>%
        mutate(data_type = data_type(ref_study)) %>%
        ggplot(aes(x = ref_study, y = auroc, fill = data_type)) +
        geom_boxplot() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        labs(x = "Reference dataset", y = "Number of reciprocal best hits",
             fill = "Data type") +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)

    # zoom in on the various classes
    classes = splitClusters(auroc, 2)
    other = classes[[1]]
    it = classes[[2]]
    glu_non_it = splitClusters(auroc[other, other], 2)[[2]]
    other = splitClusters(auroc[other, other], 2)[[1]]
    non_neuron = splitClusters(auroc[other, other], 2)[[1]]
    gaba = splitClusters(auroc[other, other], 2)[[2]]

    plotHeatmap(auroc[it, it], cex = 0.8,
                ColSideColors = dataset_to_col[getStudyId(it)],
                RowSideColors = dataset_to_col[getStudyId(it)])
    plotHeatmap(auroc[glu_non_it, glu_non_it], cex = 0.8,
                ColSideColors = dataset_to_col[getStudyId(glu_non_it)],
                RowSideColors = dataset_to_col[getStudyId(glu_non_it)])
    
    # merge L2/3 and L5 annotations
    full_data$new_labels = merge_labels(full_data$subclass_label)
    auroc = MetaNeighborUS(hvg, full_data, study_id = full_data$study_id,
                           cell_type = full_data$new_labels,
                           fast_version = TRUE)
    plotHeatmap(auroc,
                ColSideColors = dataset_to_col[getStudyId(colnames(auroc))],
                RowSideColors = dataset_to_col[getStudyId(colnames(auroc))])

    # zoom in on glutamatergic neurons
    classes = splitClusters(auroc, 3)
    it = unlist(classes[c(2,3)])
    plotHeatmap(auroc[it, it], cex = 0.8,
                ColSideColors = dataset_to_col[getStudyId(it)],
                RowSideColors = dataset_to_col[getStudyId(it)])    

    # plot # best hits and average auroc per study
    best_hits = find_best_hits(auroc) %>%
        mutate(ref_study = getStudyId(ref_cell_type))
    
    to_plot = best_hits %>%
        group_by(ref_study, ref_cell_type) %>%
        summarize(n_hits = sum(is_reciprocal)) %>%
        group_by(ref_study) %>% 
        summarize(m = mean(n_hits), sd = sd(n_hits)) %>%
        mutate(data_type = data_type(ref_study)) %>%
        ggplot(aes(x = ref_study, y = m, fill = data_type)) +
        geom_col() +
        geom_linerange(aes(ymin = m-sd, ymax = m + sd)) +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        labs(x = "Reference dataset", y = "Number of reciprocal best hits",
             fill = "Data type") +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)

    to_plot = best_hits %>%
        group_by(ref_study, ref_cell_type) %>%
        summarize(auroc = mean(auroc)) %>%
        mutate(data_type = data_type(ref_study)) %>%
        ggplot(aes(x = ref_study, y = auroc, fill = data_type)) +
        geom_boxplot() +
        theme_bw(base_size = 20) +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        labs(x = "Reference dataset", y = "Number of reciprocal best hits",
             fill = "Data type") +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)
    
    best_hits %>%
        filter(ref_study == "atac" & auroc < 0.9)

    # look at 1-vs-1 auroc
    best_hits = MetaNeighborUS(hvg, full_data, study_id = full_data$study_id,
                               cell_type = full_data$new_labels, fast_version = TRUE,
                               one_vs_best = TRUE, symmetric_output = FALSE)
    plotHeatmap(best_hits,
                ColSideColors = dataset_to_col[getStudyId(colnames(best_hits))],
                RowSideColors = dataset_to_col[getStudyId(colnames(best_hits))])
    meta_clusters = extractMetaClusters(best_hits, threshold = 0.7)
    write.table(scoreMetaClusters(meta_clusters, best_hits), "figs/multi_modal_clusters.txt")
    MetaNeighbor::scoreMetaClusters(meta_clusters, best_hits)
    print(plotUpset(meta_clusters))
    
    keep_cell = makeClusterName(full_data$study_id, full_data$subclass_label) %in% non_neuron
    auroc = MetaNeighborUS(hvg, full_data[,keep_cell],
                           study_id = full_data$study_id[keep_cell],
                           cell_type = full_data$subclass_label[keep_cell],
                           fast_version = TRUE, one_vs_best = TRUE, symmetric_output = FALSE)
    plotHeatmap(auroc, cex=0.8,
                ColSideColors = dataset_to_col[getStudyId(colnames(auroc))],
                RowSideColors = dataset_to_col[getStudyId(colnames(auroc))])
}

find_best_hits = function(auroc) {
    result = topHitsByStudy(auroc, threshold = 0, n_digits = Inf, collapse_duplicates = FALSE)
    names(result) = c("ref_cell_type", "target_cell_type", "auroc", "type")
    result = result %>%
        mutate(is_reciprocal = type == "Reciprocal_top_hit") %>%
        select(-type)
    return(result)
}

data_type = function(dataset_name) {
    result = rep("RNA", length(dataset_name))
    result[dataset_name == "atac"] = "ATAC"
    result[dataset_name == "snmc"] = "Methylation"
    return(result)
}

merge_labels = function(labels) {
    result = labels
    result[startsWith(labels, "L23")] = "L2/3 IT"
    result[startsWith(labels, "L5.IT")] = "L4/5 IT"
    result[startsWith(labels, "L5-IT")] = "L4/5 IT"
    result[startsWith(labels, "L4")] = "L4/5 IT"
    return(result)
}

if (sys.nframe() == 0) {
    main()
}
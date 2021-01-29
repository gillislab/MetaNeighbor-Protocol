
library(tidyverse)
library(SingleCellExperiment)
library(Matrix)
library(MetaNeighbor)


main = function() {
    generate_data()

    pdf("downsampling_analysis.pdf", 8, 5.5)
    make_downsampling_figures()
    dev.off()
    pdf("noise_analysis.pdf", 8, 5.5)
    make_noise_figures()
    dev.off()
}

generate_data = function() {
    # select endocrine cells
    pancreas = readRDS("../merged_pancreas.rds")
    hvg = variableGenes(pancreas, exp_labels = pancreas$study_id)
    auroc = MetaNeighborUS(hvg, pancreas, study_id = pancreas$study_id,
                           cell_type = pancreas$`cell type`, fast_version = TRUE)
    endocrine_clusters = splitClusters(auroc, 2)[[2]]
    keep_cell = makeClusterName(pancreas$study_id, pancreas$`cell type`) %in% endocrine_clusters
    endocrine = pancreas[, keep_cell]
    
    make_downsampling_analysis(endocrine, n_replicates = 10)
    make_noise_analysis(endocrine, n_replicates = 10)
}

make_downsampling_analysis = function(endocrine, n_replicates = 10) {
    # reference analysis
    static_hvg = variableGenes(endocrine, exp_labels = endocrine$study_id)
    ref_one_vs_all = MetaNeighborUS(
        static_hvg, endocrine, study_id = endocrine$study_id,
        cell_type = endocrine$`cell type`, fast_version = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(auroc_type = "1vall")
    ref_one_vs_best = MetaNeighborUS(
        static_hvg, endocrine, study_id = endocrine$study_id,
        cell_type = endocrine$`cell type`, fast_version = TRUE, one_vs_best = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(auroc_type = "1vbest")
    ref_experiments = cross_df(list(run_id = 1, f = 1, hvg = c("variable", "static"))) %>%
        full_join(bind_rows(ref_one_vs_all, ref_one_vs_best), by = character())
    
    non_baron = endocrine[, endocrine$study_id != "baron"]
    baron = endocrine[, endocrine$study_id == "baron"]
    downsample_fraction = c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001)
    experiments = cross_df(list(run_id = 1:n_replicates, f = downsample_fraction)) %>%
        group_by(run_id, f) %>%
        summarize(perform_downsampling_experiment(static_hvg, baron, non_baron, f))
    
    result = bind_rows(ref_experiments, experiments)
    write_csv(result, "downsampling_data.csv")
}

find_best_hits = function(auroc) {
    result = topHitsByStudy(auroc, threshold = 0, n_digits = Inf, collapse_duplicates = FALSE)
    names(result) = c("ref_cell_type", "target_cell_type", "auroc", "type")
    result = result %>%
        mutate(is_reciprocal = type == "Reciprocal_top_hit") %>%
        select(-type)
    return(result)
}

perform_downsampling_experiment = function(static_hvg, baron, non_baron, f) {
    new_baron = baron
    new_expr = downsample_counts(counts(baron), f)
    dimnames(new_expr) = dimnames(baron)
    counts(new_baron) = new_expr
    downsampled = cbind(new_baron, non_baron)
    hvg = variableGenes(downsampled, exp_labels = downsampled$study_id)
    if (length(hvg) < 200) {
        hvg = variableGenes(downsampled, exp_labels = downsampled$study_id, min_recurrence = 3)
    }
    
    one_vs_all = MetaNeighborUS(
        hvg, downsampled, study_id = downsampled$study_id,
        cell_type = downsampled$`cell type`, fast_version = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(hvg = "variable", auroc_type = "1vall")
    
    static_one_vs_all = MetaNeighborUS(
        static_hvg, downsampled, study_id = downsampled$study_id,
        cell_type = downsampled$`cell type`, fast_version = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(hvg = "static", auroc_type = "1vall")
    
    one_vs_best = MetaNeighborUS(
        hvg, downsampled, study_id = downsampled$study_id,
        cell_type = downsampled$`cell type`, fast_version = TRUE, one_vs_best = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(hvg = "variable", auroc_type = "1vbest")
    
    static_one_vs_best = MetaNeighborUS(
        static_hvg, downsampled, study_id = downsampled$study_id,
        cell_type = downsampled$`cell type`, fast_version = TRUE, one_vs_best = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(hvg = "static", auroc_type = "1vbest")
    
    result = bind_rows(one_vs_all, static_one_vs_all,
                       one_vs_best, static_one_vs_best)
    return(result)
}

downsample_counts = function(expr_matrix, f = 0.1) {
    result = apply(expr_matrix, 2, downsample_cell, f)
    return(Matrix(result, sparse = TRUE))
}

downsample_cell = function(expr, f = 0.1) {
    total_counts = ceiling(sum(expr))
    result = rmultinom(n = 1, size = f*total_counts, prob = expr/total_counts)
    return(result)
}

make_downsampling_figures = function() {
    best_hits = read_csv("downsampling_data.csv")
    plot_downsampling_analysis(filter(best_hits, auroc_type == "1vall"))
    plot_downsampling_analysis(filter(best_hits, auroc_type == "1vbest"), 0.7)
}

plot_downsampling_analysis = function(best_hits, auroc_threshold = 0.9, mean_umis = 5000) {
    best_hits = best_hits %>%
        mutate(downsampling_number = as.factor(as.numeric(f) * mean_umis)) %>%
        mutate(downsampling_number = fct_rev(downsampling_number))

    to_plot = best_hits %>%
        group_by(downsampling_number, ref_cell_type, hvg, run_id) %>%
        summarize(mean_auroc = mean(auroc))  %>%
        group_by(downsampling_number, ref_cell_type, hvg) %>%
        summarize(mean_auroc = mean(mean_auroc))  %>%
        ggplot(aes(x = downsampling_number, y = mean_auroc,
                   col = ref_cell_type, linetype = hvg)) +
        geom_point() +
        geom_line(aes(group = paste(ref_cell_type, hvg))) +
        geom_hline(yintercept = auroc_threshold, linetype = "dashed") +
        theme_bw(base_size = 20) +
        labs(x = "Median UMIs per cell after downsampling",
             y = "Average AUROC",
             col = "Cell type",
             linetype = "HVG selection") +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        scale_color_brewer(palette = "Set1")
    print(to_plot)
    
    to_plot = best_hits %>%
        filter(hvg == "variable") %>%
        group_by(downsampling_number, ref_cell_type, run_id) %>%
        summarize(n_reciprocal_hits = sum(is_reciprocal)) %>%
        group_by(downsampling_number, ref_cell_type) %>%
        summarize(n_reciprocal_hits = mean(n_reciprocal_hits)) %>%
        ggplot(aes(x = downsampling_number, y = n_reciprocal_hits, fill = ref_cell_type)) +
        geom_col() +
        theme_bw(base_size = 20) +
        labs(x = "Median UMIs per cell after downsampling",
             y= "Number of reciprocal top hits",
             fill = "Cell type") +
        theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)
}


make_noise_analysis = function(endocrine, n_replicates = 10) {
    hvg = variableGenes(endocrine, exp_labels = endocrine$study_id)
    non_baron = endocrine[hvg, endocrine$study_id != "baron"]
    baron = endocrine[hvg, endocrine$study_id == "baron"]
    counts(non_baron) = as.matrix(counts(non_baron))
    counts(baron) = as.matrix(counts(baron))
    
    noise_level = seq(0, 5, by = 0.25)
    result = cross_df(list(run_id = 1:n_replicates, d = noise_level)) %>%
        group_by(run_id, d) %>%
        summarize(perform_noise_experiment(hvg, baron, non_baron, d))
    write_csv(result, "noise_data.csv")
}

perform_noise_experiment = function(hvg, baron, non_baron, d) {
    new_baron = baron
    new_expr = make_artificial_batch(counts(baron), d)
    dimnames(new_expr) = dimnames(baron)
    counts(new_baron) = new_expr
    noisy_data = cbind(new_baron, non_baron)

    one_vs_all = MetaNeighborUS(
        hvg, noisy_data, study_id = noisy_data$study_id,
        cell_type = noisy_data$`cell type`, fast_version = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(auroc_type = "1vall")

    one_vs_best = MetaNeighborUS(
        hvg, noisy_data, study_id = noisy_data$study_id,
        cell_type = noisy_data$`cell type`, fast_version = TRUE, one_vs_best = TRUE
    ) %>%
        find_best_hits() %>%
        filter(getStudyId(ref_cell_type) == "baron") %>%
        add_column(auroc_type = "1vbest")

    result = bind_rows(one_vs_all, one_vs_best)
    return(result)
}

# see https://www.biorxiv.org/content/10.1101/2020.05.22.111211v2.full, "Artifical batches"
make_artificial_batch = function(counts, d) {
    counts = as.matrix(counts)
    noise = rnorm(length(counts), mean = 0, sd = d)
    result = counts + noise
    result[result<0] = 0
    return(result)
}

make_noise_figures = function() {
    best_hits = read_csv("noise_data.csv") %>%
        dplyr::rename(noise_level = d)
    plot_noise_analysis(filter(best_hits, auroc_type == "1vall"))
    plot_noise_analysis(filter(best_hits, auroc_type == "1vbest"), 0.7)
}

plot_noise_analysis = function(best_hits, auroc_threshold = 0.9, mean_umi = 17) {
    to_plot = best_hits %>%
        group_by(noise_level, ref_cell_type, run_id) %>%
        summarize(mean_auroc = mean(auroc))  %>%
        group_by(noise_level, ref_cell_type) %>%
        summarize(mean_auroc = mean(mean_auroc))  %>%
        ggplot(aes(x = noise_level / mean_umi, y = mean_auroc, col = ref_cell_type)) +
        geom_point() +
        geom_line() +
        theme_bw(base_size = 20) +
        labs(x = "Batch effect size (noise level)", y = "Average AUROC",
             col = "Cell type") +
        geom_hline(yintercept = auroc_threshold, linetype = "dashed") +
        scale_x_continuous(labels = scales::percent) +
        scale_color_brewer(palette = "Set1")
    print(to_plot)
    
    to_plot = best_hits %>%
        group_by(noise_level, ref_cell_type, run_id) %>%
        summarize(n_reciprocal_hits = sum(is_reciprocal)) %>%
        group_by(noise_level, ref_cell_type) %>%
        summarize(n_reciprocal_hits = mean(n_reciprocal_hits)) %>%
        ggplot(aes(x = noise_level / mean_umi, y = n_reciprocal_hits, fill = ref_cell_type)) +
        geom_col() +
        theme_bw(base_size = 20) +
        labs(x = "Batch effect size (noise level)",
             y= "Number of reciprocal top hits",
             fill = "Cell type") +
        scale_x_continuous(labels = scales::percent) +
        scale_fill_brewer(palette = "Set1")
    print(to_plot)
}

if (sys.nframe() == 0) {
    main()
}
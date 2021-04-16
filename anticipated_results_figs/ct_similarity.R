# equivalent to code shown in anticipated results section, but with batch annotations shown.

library(MetaNeighbor)
library(SingleCellExperiment)

main = function() {
    biccn_data = readRDS("../data/biccn_gaba.rds")
    biccn_hvgs = variableGenes(biccn_data, exp_labels = biccn_data$study_id)

    cell_types = as.factor(
        makeClusterName(biccn_data$study_id, biccn_data$joint_subclass_label)
    )
    normalization_factor = Matrix::colSums(assay(biccn_data)) / 1000000
    cpm = assay(biccn_data)
    cpm@x = cpm@x / rep.int(normalization_factor, diff(cpm@p))
    cpm = as.matrix(cpm[biccn_hvgs,])

    centroids = sapply(levels(cell_types), function(ct) {        
        matrixStats::rowMeans2(log2(cpm+1), cols = cell_types == ct)
    })

    centroid_cor = cor(centroids, method = "spearman")
    aurocs = MetaNeighborUS(var_genes = biccn_hvgs,
                            dat = biccn_data,
                            study_id = biccn_data$study_id,
                            cell_type = biccn_data$joint_subclass_label,
                            fast_version = TRUE)
    
    dataset_cols = c(scCv2="#2A2F7A", scCv3="#4E517A",
                     snCv2="#A0BCFB", snCv3M="#5586FA", snCv3Z="#4B77DE",
                     scSS="#C73733", snSS="#FB908C")
    ct_names = sort(unique(biccn_data$joint_subclass_label))
    ct_cols = RColorBrewer::brewer.pal(n = length(ct_names), "Set1")
    names(ct_cols) = ct_names
    
    pdf("ct_similarity.pdf")
    plotHeatmap((1+centroid_cor)/2, cex=0.9, labCol=FALSE, dendrogram="col",
                ColSideColors = dataset_cols[getStudyId(colnames(centroid_cor))],
                RowSideColors = ct_cols[getCellType(colnames(centroid_cor))])
    legend("topleft", legend = names(dataset_cols), fill = dataset_cols, inset = c(-0.1,0.2),
           title = "Dataset", cex = 0.6)
    legend("topleft", legend = names(ct_cols), fill = ct_cols, inset = c(-0.1,0.5),
           title = "Cell type", cex = 0.6)
    
    plotHeatmap(aurocs, cex=0.9, labCol=FALSE, dendrogram="col",
                ColSideColors = dataset_cols[getStudyId(colnames(aurocs))],
                RowSideColors = ct_cols[getCellType(colnames(aurocs))])
    legend("topleft", legend = names(dataset_cols), fill = dataset_cols, inset = c(-0.1,0.2),
           title = "Dataset", cex = 0.6)
    legend("topleft", legend = names(ct_cols), fill = ct_cols, inset = c(-0.1,0.5),
           title = "Cell type", cex = 0.6)
    dev.off()
}

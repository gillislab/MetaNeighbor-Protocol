
library(tidyverse)


base_url = "http://data.nemoarchive.org/biccn/lab/ecker/chromatin/scell/processed/analysis/EckerRen_Mouse_MOp_methylation_ATAC"
snatac_metadata = file.path(base_url, "study/ATAC/MOp.snATAC-seq.AnalysisResult.csv.gz")
snatac_data_dir = file.path(base_url, "dataset/ATAC")
snmc_metadata = file.path(base_url, "study/mC/MOp_clustering/MOp.snmC-seq.AnalysisResult.csv.gz")
snmc_data_dir = file.path(base_url, "dataset/mC")

# Download snatac
output_dir = "snatac"
dir.create(output_dir, showWarnings = FALSE)
metadata_file = file.path(output_dir, "metadata.csv.gz")
download.file(snatac_metadata, metadata_file)
metadata = read_csv(metadata_file)
sample_names = unique(metadata$sample)
for (sample_name in sample_names) {
    filename = paste0(sample_name, ".snap")
    download.file(file.path(snatac_data_dir, filename),
                  file.path(output_dir, filename))
}

# Download snmC
output_dir = "snmc"
dir.create(output_dir, showWarnings = FALSE)
metadata_file = file.path(output_dir, "metadata.csv.gz")
download.file(snmc_metadata, metadata_file)
metadata = read_csv(metadata_file)
sample_names = unique(paste(metadata$Region, metadata$FACS_Date, sep="-"))
for (sample_name in sample_names) {
    filename = paste0(sample_name, ".mcds")
    download.file(file.path(snmc_data_dir, filename),
                  file.path(output_dir, filename))
}

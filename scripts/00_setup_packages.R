# =============================================================
# 00. Package setup (run once per computer)
# =============================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("dada2", "ShortRead", "Biostrings", "phyloseq",
                       "DECIPHER", "ANCOMBC"), update = FALSE, ask = FALSE)
install.packages(c("phangorn", "picante", "vegan", "ggplot2",
                   "RColorBrewer", "patchwork"))

# Raw data: run scripts/02_download_ena.R. Two runs need the NCBI SRA Toolkit
# (fasterq-dump); see the notes at the top of that script.
# =============================================================
# 04. Quality filtering and trimming of paired-end reads
# Primers were already removed by the submitters (Episode 03).
# V4 amplicon (515F-806R) is ~253 bp; reads are 247 bp.
# =============================================================

suppressPackageStartupMessages({
  library(dada2)
  library(ggplot2)
})

meta <- read.csv("data/sample-metadata.csv")
fwd_files <- file.path("data/raw", paste0(meta$run, "_1.fastq.gz"))
rev_files <- file.path("data/raw", paste0(meta$run, "_2.fastq.gz"))

amplicon_len <- 253

# ---- 1. Quick test of truncation settings (4 samples) ----
test_idx <- match(unique(meta$compartment), meta$compartment)   # one per compartment
settings <- list(c(240, 180), c(230, 160), c(220, 140))
test_dir <- file.path(tempdir(), "trunc_test")
dir.create(test_dir, showWarnings = FALSE)

trunc_test <- do.call(rbind, lapply(settings, function(tl) {
  out <- filterAndTrim(fwd_files[test_idx], file.path(test_dir, paste0("F_", test_idx, ".fq.gz")),
                       rev_files[test_idx], file.path(test_dir, paste0("R_", test_idx, ".fq.gz")),
                       truncLen = tl, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                       rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  data.frame(truncLen = paste(tl, collapse = "/"),
             overlap_bp = sum(tl) - amplicon_len,
             pct_kept = round(100 * sum(out[, "reads.out"]) / sum(out[, "reads.in"]), 1))
}))
cat("Truncation test (4 samples):\n"); print(trunc_test, row.names = FALSE)
unlink(test_dir, recursive = TRUE)

# ---- 2. Filter all samples with the chosen setting ----
trunc_len <- c(230, 160)

filt_dir <- "data/processed/filtered"
dir.create(filt_dir, recursive = TRUE, showWarnings = FALSE)
filt_f <- file.path(filt_dir, paste0(meta$run, "_F_filt.fastq.gz"))
filt_r <- file.path(filt_dir, paste0(meta$run, "_R_filt.fastq.gz"))

out <- filterAndTrim(fwd_files, filt_f, rev_files, filt_r,
                     truncLen = trunc_len, maxN = 0, maxEE = c(2, 2), truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, multithread = FALSE)

# ---- 3. Summary ----
meta$reads_in  <- out[, "reads.in"]
meta$reads_out <- out[, "reads.out"]
meta$pct_kept  <- round(100 * meta$reads_out / meta$reads_in, 1)

cat("\nTotal pairs in: ", format(sum(meta$reads_in),  big.mark = ","),
    "\nTotal pairs out:", format(sum(meta$reads_out), big.mark = ","),
    sprintf("\nOverall kept:    %.1f%%\n", 100 * sum(meta$reads_out) / sum(meta$reads_in)))
print(summary(meta$pct_kept))
cat("\nMean % kept by compartment:\n")
print(aggregate(pct_kept ~ compartment, data = meta, FUN = function(x) round(mean(x), 1)))
cat("\nLowest-retention samples:\n")
print(head(meta[order(meta$pct_kept), c("sample_id", "compartment", "day",
                                        "reads_in", "reads_out", "pct_kept")], 5),
      row.names = FALSE)

write.csv(meta[, c("run", "sample_id", "compartment", "day",
                   "reads_in", "reads_out", "pct_kept")],
          "results/tables/04_filter_summary.csv", row.names = FALSE)

# ---- 4. Quality after filtering (one sample per compartment) ----
p_ff <- plotQualityProfile(filt_f[test_idx]) + ggtitle("Forward reads after filtering")
p_rf <- plotQualityProfile(filt_r[test_idx]) + ggtitle("Reverse reads after filtering")
print(p_ff); print(p_rf)
ggsave("results/figures/04_quality_forward_filtered.png", p_ff, width = 10, height = 7, dpi = 300)
ggsave("results/figures/04_quality_reverse_filtered.png", p_rf, width = 10, height = 7, dpi = 300)
# =============================================================
# 05. DADA2: error models, denoising, merging pairs, chimeras
# Resumable: each sample is saved as it finishes; re-Source to continue.
# =============================================================

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
  library(ggplot2)
})
set.seed(100)

meta   <- read.csv("data/sample-metadata.csv")
filt_f <- file.path("data/processed/filtered", paste0(meta$run, "_F_filt.fastq.gz"))
filt_r <- file.path("data/processed/filtered", paste0(meta$run, "_R_filt.fastq.gz"))
names(filt_f) <- names(filt_r) <- meta$sample_id

# ---- 1. Error models (saved, so learned only once) ----
err_file <- "data/processed/error_models.rds"
if (file.exists(err_file)) {
  errs <- readRDS(err_file)
} else {
  errs <- list(F = learnErrors(filt_f, multithread = FALSE),
               R = learnErrors(filt_r, multithread = FALSE))
  saveRDS(errs, err_file)
}
p_errF <- plotErrors(errs$F, nominalQ = TRUE) + ggtitle("Error model: forward reads")
p_errR <- plotErrors(errs$R, nominalQ = TRUE) + ggtitle("Error model: reverse reads")
print(p_errF); print(p_errR)
ggsave("results/figures/05_error_model_forward.png", p_errF, width = 10, height = 8, dpi = 300)
ggsave("results/figures/05_error_model_reverse.png", p_errR, width = 10, height = 8, dpi = 300)

# ---- 2. Denoise and merge, one sample at a time (resumable) ----
ckpt_dir <- "data/processed/dada_per_sample"
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

for (s in meta$sample_id) {
  ckpt <- file.path(ckpt_dir, paste0(s, ".rds"))
  if (file.exists(ckpt)) next
  t0 <- Sys.time()
  dF <- dada(filt_f[[s]], err = errs$F, multithread = FALSE, verbose = FALSE)
  dR <- dada(filt_r[[s]], err = errs$R, multithread = FALSE, verbose = FALSE)
  mg <- mergePairs(dF, filt_f[[s]], dR, filt_r[[s]], minOverlap = 20)
  saveRDS(list(denoisedF = sum(getUniques(dF)),
               denoisedR = sum(getUniques(dR)),
               merged    = mg), ckpt)
  cat(sprintf("%s done (%d/%d) in %.1f min\n", s, match(s, meta$sample_id),
              nrow(meta), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

res <- lapply(meta$sample_id, function(s) readRDS(file.path(ckpt_dir, paste0(s, ".rds"))))
names(res) <- meta$sample_id

# ---- 3. Sequence table and length check ----
seqtab <- makeSequenceTable(lapply(res, `[[`, "merged"))
cat("\nSequence table (samples x ASVs):", dim(seqtab), "\n")

len_tab <- tapply(colSums(seqtab), nchar(colnames(seqtab)), sum)
len_pct <- round(100 * len_tab / sum(len_tab), 2)
cat("\nMerged length (bp): % of reads (lengths with >= 0.1% shown)\n")
print(len_pct[len_pct >= 0.1])

keep_len <- 250:256   # expected V4 (515F-806R) length range
seqtab_len <- seqtab[, nchar(colnames(seqtab)) %in% keep_len]
cat(sprintf("\nReads within %d-%d bp: %.1f%%\n", min(keep_len), max(keep_len),
            100 * sum(seqtab_len) / sum(seqtab)))

# ---- 4. Remove chimeras ----
seqtab_nochim <- removeBimeraDenovo(seqtab_len, method = "consensus",
                                    multithread = FALSE, verbose = TRUE)
cat("ASVs after chimera removal:", ncol(seqtab_nochim), "\n")
cat(sprintf("Reads kept after chimera removal: %.1f%%\n",
            100 * sum(seqtab_nochim) / sum(seqtab_len)))

# ---- 5. Track reads through the pipeline ----
filt_sum <- read.csv("results/tables/04_filter_summary.csv")
track <- data.frame(
  sample_id    = meta$sample_id,
  compartment  = meta$compartment,
  day          = meta$day,
  raw          = filt_sum$reads_in[match(meta$sample_id, filt_sum$sample_id)],
  filtered     = filt_sum$reads_out[match(meta$sample_id, filt_sum$sample_id)],
  denoised_F   = sapply(res, `[[`, "denoisedF"),
  denoised_R   = sapply(res, `[[`, "denoisedR"),
  merged       = sapply(res, function(x) sum(getUniques(x$merged))),
  length_ok    = rowSums(seqtab_len)[meta$sample_id],
  non_chimeric = rowSums(seqtab_nochim)[meta$sample_id]
)
track$pct_raw_kept <- round(100 * track$non_chimeric / track$raw, 1)
track$pct_merged   <- round(100 * track$merged / track$filtered, 1)

cat("\nRead tracking summary:\n")
print(summary(track[, c("non_chimeric", "pct_raw_kept", "pct_merged")]))
cat("\nMean % of raw pairs kept, by compartment:\n")
print(aggregate(pct_raw_kept ~ compartment, data = track, FUN = function(x) round(mean(x), 1)))
cat("\nLowest final depth:\n")
print(head(track[order(track$non_chimeric),
                 c("sample_id", "compartment", "day", "merged", "non_chimeric", "pct_raw_kept")], 5),
      row.names = FALSE)
write.csv(track, "results/tables/05_read_tracking.csv", row.names = FALSE)

# ---- 6. Save outputs ----
asv_seqs <- colnames(seqtab_nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
writeXStringSet(DNAStringSet(setNames(asv_seqs, asv_ids)), "results/05_asv_sequences.fasta")
saveRDS(seqtab_nochim, "data/processed/seqtab_nochim.rds")
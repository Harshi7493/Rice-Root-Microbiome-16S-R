# =============================================================
# 06. Taxonomic classification (SILVA 138.2, DADA2 naive Bayes)
# Database: Callahan B. (2024) Silva taxonomic training data formatted
#   for DADA2 (v138.2). Zenodo. doi:10.5281/zenodo.14169026
# ASVs are classified in chunks of 500 to limit memory use
# (results are identical); finished chunks are saved, so the script
# can be re-Sourced to resume after an interruption.
# =============================================================

suppressPackageStartupMessages(library(dada2))
options(timeout = 1800)

# ---- 1. SILVA reference files (copy from first project, else download) ----
ref_dir <- "data/reference"
dir.create(ref_dir, recursive = TRUE, showWarnings = FALSE)
refs <- c(genus   = "silva_nr99_v138.2_toGenus_trainset.fa.gz",
          species = "silva_v138.2_assignSpecies.fa.gz")
md5_expected <- c(genus   = "1764e2a36b4500ccb1c7d5261948a414",
                  species = "62fd939bd1aa9832d2089847f9f0b87d")
old_ref <- "C:/Users/Harshi/Documents/16S-Amplicon-Microbiome-Analysis-R/data/reference"

for (r in names(refs)) {
  dest <- file.path(ref_dir, refs[[r]])
  if (!file.exists(dest) && file.exists(file.path(old_ref, refs[[r]]))) {
    file.copy(file.path(old_ref, refs[[r]]), dest)
  }
  if (!file.exists(dest)) {
    download.file(paste0("https://zenodo.org/records/14169026/files/", refs[[r]],
                         "?download=1"), dest, mode = "wb")
  }
  ok <- unname(tools::md5sum(dest)) == md5_expected[[r]]
  cat(refs[[r]], "- checksum OK:", ok, "\n")
  if (!ok) stop("Checksum mismatch: delete ", dest, " and run again.")
}

# ---- 2. Assign taxonomy in chunks (resumable) ----
seqtab_nochim <- readRDS("data/processed/seqtab_nochim.rds")
asv_seqs <- colnames(seqtab_nochim)
asv_ids  <- paste0("ASV", seq_along(asv_seqs))
cat("ASVs to classify:", length(asv_seqs), "\n")

chunk_size <- 500
chunks <- split(seq_along(asv_seqs), ceiling(seq_along(asv_seqs) / chunk_size))
ckpt_dir <- "data/processed/taxonomy_chunks"
dir.create(ckpt_dir, recursive = TRUE, showWarnings = FALSE)

for (k in seq_along(chunks)) {
  ckpt <- file.path(ckpt_dir, sprintf("chunk_%02d.rds", k))
  if (file.exists(ckpt)) next
  t0 <- Sys.time()
  set.seed(100 + k)   # reproducible bootstrap per chunk
  tx <- assignTaxonomy(asv_seqs[chunks[[k]]], file.path(ref_dir, refs[["genus"]]),
                       minBoot = 50, multithread = FALSE)
  tx <- addSpecies(tx, file.path(ref_dir, refs[["species"]]))
  saveRDS(tx, ckpt)
  rm(tx); invisible(gc())   # free memory before the next chunk
  cat(sprintf("Chunk %d/%d done in %.1f min\n", k, length(chunks),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

taxa <- do.call(rbind, lapply(seq_along(chunks), function(k)
  readRDS(file.path(ckpt_dir, sprintf("chunk_%02d.rds", k)))))
stopifnot(identical(rownames(taxa), asv_seqs))   # same order as the ASV table

# ---- 3. Summaries (by ASVs and by reads) ----
reads_per_asv <- colSums(seqtab_nochim)

cat("\n% of ASVs classified at each rank:\n")
print(round(100 * colMeans(!is.na(taxa)), 1))

cat("\nKingdoms (ASVs and % of reads):\n")
k <- ifelse(is.na(taxa[, "Kingdom"]), "Unassigned", taxa[, "Kingdom"])
print(data.frame(ASVs = as.vector(table(k)),
                 pct_reads = round(100 * tapply(reads_per_asv, k, sum) / sum(reads_per_asv), 2),
                 row.names = names(table(k))))

is_chloro <- taxa[, "Order"]  %in% "Chloroplast"
is_mito   <- taxa[, "Family"] %in% "Mitochondria"
cat(sprintf("\nChloroplast:  %d ASVs, %.1f%% of reads\n",
            sum(is_chloro), 100 * sum(reads_per_asv[is_chloro]) / sum(reads_per_asv)))
cat(sprintf("Mitochondria: %d ASVs, %.1f%% of reads\n",
            sum(is_mito), 100 * sum(reads_per_asv[is_mito]) / sum(reads_per_asv)))

# What makes up the 251-bp group seen in Episode 05?
len <- nchar(asv_seqs)
grp <- ifelse(is_chloro, "Chloroplast", ifelse(is_mito, "Mitochondria",
                                               ifelse(is.na(taxa[, "Kingdom"]), "Unassigned", taxa[, "Kingdom"])))
cat("\n% of reads by length and group:\n")
print(round(100 * tapply(reads_per_asv, list(len, grp), sum) / sum(reads_per_asv), 2))

cat("\nTop 10 phyla by % of reads:\n")
ph <- ifelse(is.na(taxa[, "Phylum"]), "Unassigned", taxa[, "Phylum"])
print(head(sort(round(100 * tapply(reads_per_asv, ph, sum) / sum(reads_per_asv), 1),
                decreasing = TRUE), 10))

# Archaea, including methanogens
arch <- taxa[, "Kingdom"] %in% "Archaea"
cat("\nArchaeal genera (ASVs):\n")
print(head(sort(table(taxa[arch, "Genus"], useNA = "ifany"), decreasing = TRUE), 8))

# ---- 4. Save ----
taxa_table <- data.frame(ASV = asv_ids, taxa, length = len, total_reads = reads_per_asv,
                         Sequence = asv_seqs, row.names = NULL)
write.csv(taxa_table, "results/tables/06_taxonomy.csv", row.names = FALSE)
saveRDS(taxa, "data/processed/taxonomy.rds")
# =============================================================
# 02. Download raw paired-end FASTQ files
# Study: Edwards et al. (2015) PNAS 112:E911-E920, SRA SRP044745
# Subset: time-series experiment, M104 in Davis soil,
#         days 1, 5, 13 x 4 compartments x 3 replicates (36 runs)
#
# 34 runs are downloaded from ENA and verified by MD5 checksum.
# 2 runs (SRR2107429, SRR2107479) are missing from ENA's file server,
# so they were downloaded from NCBI SRA with the SRA Toolkit
# (fasterq-dump --split-files) and verified by matching
# forward/reverse read counts (116,739 and 93,373 read pairs).
# =============================================================

options(timeout = 3600)   # allow time for large files

meta <- read.csv("data/sample-metadata.csv")
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# ---- 1. Ask ENA for file links, sizes and checksums ----
ena_query <- function(run) {
  url <- paste0("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=", run,
                "&result=read_run&fields=run_accession,library_layout,",
                "fastq_ftp,fastq_md5,fastq_bytes&format=tsv")
  read.delim(url, colClasses = "character")
}

manifest <- do.call(rbind, lapply(meta$run, function(r) {
  x <- ena_query(r)
  data.frame(run    = r,
             layout = x$library_layout,
             url    = strsplit(x$fastq_ftp, ";")[[1]],
             md5    = strsplit(x$fastq_md5, ";")[[1]],
             bytes  = as.numeric(strsplit(x$fastq_bytes, ";")[[1]]))
}))
manifest$file <- basename(manifest$url)

cat("Files to download:", nrow(manifest),
    sprintf("(%.0f MB in total)\n", sum(manifest$bytes) / 1e6))
print(table(files_per_run = table(manifest$run)))

# ---- 2. Download and verify (skips files already downloaded) ----
ncbi_runs <- c("SRR2107429", "SRR2107479")   # fetched from NCBI SRA instead

for (i in seq_len(nrow(manifest))) {
  dest <- file.path("data/raw", manifest$file[i])
  from_ncbi <- manifest$run[i] %in% ncbi_runs
  is_ok <- function() file.exists(dest) &&
    (from_ncbi || unname(tools::md5sum(dest)) == manifest$md5[i])
  
  ok <- is_ok()
  attempt <- 0
  while (!ok && !from_ncbi && attempt < 3) {
    attempt <- attempt + 1
    cat(sprintf("[%d/%d] %s (attempt %d)\n", i, nrow(manifest),
                manifest$file[i], attempt))
    try(download.file(paste0("https://", manifest$url[i]), dest,
                      mode = "wb", quiet = TRUE))
    ok <- is_ok()
  }
  manifest$md5_ok[i] <- ok
  manifest$source[i] <- if (from_ncbi) "NCBI SRA (fasterq-dump)" else "ENA"
}

# ---- 3. Summary ----
cat("\nAll checksums OK:", all(manifest$md5_ok), "\n")
if (!all(manifest$md5_ok)) print(manifest[!manifest$md5_ok, c("run", "file")])
print(table(manifest$source) / 2)   # number of runs from each source

write.csv(manifest, "data/ena-manifest.csv", row.names = FALSE)
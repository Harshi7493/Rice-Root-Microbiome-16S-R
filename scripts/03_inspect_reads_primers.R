# =============================================================
# 03. Inspect raw paired-end reads and detect PCR primers
# =============================================================

suppressPackageStartupMessages({
  library(ShortRead)
  library(Biostrings)
  library(dada2)
  library(ggplot2)
})
set.seed(100)

meta <- read.csv("data/sample-metadata.csv")
fwd_files <- file.path("data/raw", paste0(meta$run, "_1.fastq.gz"))
rev_files <- file.path("data/raw", paste0(meta$run, "_2.fastq.gz"))
stopifnot(all(file.exists(fwd_files)), all(file.exists(rev_files)))

# ---- 1. Reads per sample ----
meta$reads_fwd <- sapply(fwd_files, function(f) countFastq(f)$records)
meta$reads_rev <- sapply(rev_files, function(f) countFastq(f)$records)
cat("Forward and reverse counts identical for every sample:",
    all(meta$reads_fwd == meta$reads_rev), "\n")
cat("Total read pairs:", format(sum(meta$reads_fwd), big.mark = ","), "\n")
print(summary(meta$reads_fwd))
print(aggregate(reads_fwd ~ compartment, data = meta, FUN = median))

# ---- 2. Read lengths (sample of 5 files) ----
idx <- sample(seq_along(fwd_files), 5)
len_f <- unlist(lapply(fwd_files[idx], function(f) width(readFastq(f))))
len_r <- unlist(lapply(rev_files[idx], function(f) width(readFastq(f))))
cat("\nForward read lengths:\n"); print(head(sort(table(len_f), decreasing = TRUE), 3))
cat("Reverse read lengths:\n");   print(head(sort(table(len_r), decreasing = TRUE), 3))

# ---- 3. Detect primers ----
# Candidate 16S primers for V4 / V4-V5 / V3-V4 regions
primers <- c(
  "515F (GTGCCAGCMGCCGCGGTAA)"      = "GTGCCAGCMGCCGCGGTAA",
  "515F-Y (GTGYCAGCMGCCGCGGTAA)"    = "GTGYCAGCMGCCGCGGTAA",
  "341F (CCTACGGGNGGCWGCAG)"        = "CCTACGGGNGGCWGCAG",
  "806R (GGACTACHVGGGTWTCTAAT)"     = "GGACTACHVGGGTWTCTAAT",
  "907R (CCGTCAATTCMTTTRAGTTT)"     = "CCGTCAATTCMTTTRAGTTT",
  "926R (CCGYCAATTYMTTTRAGTTT)"     = "CCGYCAATTYMTTTRAGTTT",
  "785R (GACTACHVGGGTATCTAATCC)"    = "GACTACHVGGGTATCTAATCC"
)

# 20,000 reads from 6 samples spread across compartments
pick <- unlist(lapply(split(seq_len(nrow(meta)), meta$compartment), `[`, 1:2))[1:6]
sample_reads <- function(files) {
  do.call(c, lapply(files[pick], function(f) {
    s <- FastqSampler(f, n = 20000 / length(pick)); x <- yield(s); close(s); sread(x)
  }))
}
reads_f <- sample_reads(fwd_files)
reads_r <- sample_reads(rev_files)

primer_scan <- function(reads, read_label) {
  heads <- subseq(reads, 1, pmin(width(reads), 60))   # look in the first 60 bases
  do.call(rbind, lapply(names(primers), function(p) {
    m <- vmatchPattern(DNAString(primers[[p]]), heads,
                       max.mismatch = 1, fixed = FALSE)
    starts <- sapply(startIndex(m), function(s) if (length(s)) s[1] else NA)
    data.frame(read = read_label, primer = p,
               pct_reads_with_primer = round(100 * mean(!is.na(starts)), 1),
               most_common_start = if (all(is.na(starts))) NA
               else as.integer(names(which.max(table(starts)))))
  }))
}
primer_hits <- rbind(primer_scan(reads_f, "forward (R1)"),
                     primer_scan(reads_r, "reverse (R2)"))
cat("\nPrimer detection (first 60 bases of each read):\n")
print(primer_hits[order(primer_hits$read, -primer_hits$pct_reads_with_primer), ],
      row.names = FALSE)

# Most common first 25 bases, to see what precedes / forms the primer
cat("\nMost common read starts (first 25 bases):\n")
cat("Forward:\n"); print(head(sort(table(as.character(subseq(reads_f, 1, 25))),
                                   decreasing = TRUE), 3))
cat("Reverse:\n"); print(head(sort(table(as.character(subseq(reads_r, 1, 25))),
                                   decreasing = TRUE), 3))

write.csv(meta, "results/tables/03_reads_per_sample.csv", row.names = FALSE)
write.csv(primer_hits, "results/tables/03_primer_detection.csv", row.names = FALSE)

# ---- 4. Quality profiles (one sample per compartment) ----
one_per_comp <- match(unique(meta$compartment), meta$compartment)
p_f <- plotQualityProfile(fwd_files[one_per_comp]) + ggtitle("Forward reads (R1)")
p_r <- plotQualityProfile(rev_files[one_per_comp]) + ggtitle("Reverse reads (R2)")
print(p_f); print(p_r)
ggsave("results/figures/03_quality_forward.png", p_f, width = 10, height = 7, dpi = 300)
ggsave("results/figures/03_quality_reverse.png", p_r, width = 10, height = 7, dpi = 300)
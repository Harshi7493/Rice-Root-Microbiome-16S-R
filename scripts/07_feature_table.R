# =============================================================
# 07. Phyloseq object, organelle analysis and non-target removal
# Reproduces the microbial:organelle read analysis of
# Edwards et al. (2015) Fig. 6A, then removes host sequences.
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(Biostrings)
  library(ggplot2)
})

# ---- 1. Build the phyloseq object ----
seqtab <- readRDS("data/processed/seqtab_nochim.rds")
taxa   <- readRDS("data/processed/taxonomy.rds")
meta   <- read.csv("data/sample-metadata.csv")
rownames(meta) <- meta$sample_id
meta$compartment <- factor(meta$compartment,
                           levels = c("Bulk soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
meta$day_f <- factor(paste("Day", meta$day), levels = c("Day 1", "Day 5", "Day 13"))
stopifnot(all(rownames(seqtab) %in% rownames(meta)))

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(meta))
dna <- DNAStringSet(taxa_names(ps)); names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq_len(ntaxa(ps)))
ps

# ---- 2. Classify each ASV ----
tax <- as.data.frame(tax_table(ps))
category <- ifelse(tax$Family %in% "Mitochondria", "Mitochondria",
                   ifelse(tax$Order  %in% "Chloroplast",  "Chloroplast",
                          ifelse(tax$Kingdom %in% c("Bacteria", "Archaea") & !is.na(tax$Phylum),
                                 "Microbial", "Other non-target")))
cat("ASVs per category:\n"); print(table(category))

# ---- 3. Organelle share per sample (cf. Edwards et al. Fig. 6A) ----
otu <- as(otu_table(ps), "matrix")
cat_reads <- sapply(c("Microbial", "Chloroplast", "Mitochondria", "Other non-target"),
                    function(cg) rowSums(otu[, category == cg, drop = FALSE]))
organelle <- data.frame(meta[rownames(otu), c("sample_id", "compartment", "day", "day_f")],
                        round(100 * cat_reads / rowSums(otu), 2), check.names = FALSE)
organelle$organelle_pct <- organelle$Chloroplast + organelle$Mitochondria

cat("\nMean % of reads by category and compartment:\n")
print(aggregate(cbind(Microbial, Chloroplast, Mitochondria) ~ compartment,
                data = organelle, FUN = function(x) round(mean(x), 1)))
cat("\nEndosphere: mean % microbial reads by day:\n")
print(aggregate(Microbial ~ day_f, data = organelle[organelle$compartment == "Endosphere", ],
                FUN = function(x) round(mean(x), 1)))
write.csv(organelle, "results/tables/07_organelle_by_sample.csv", row.names = FALSE)

p_org <- ggplot(organelle, aes(day_f, Microbial, colour = compartment)) +
  geom_point(size = 3, alpha = 0.8, position = position_dodge(width = 0.4)) +
  stat_summary(aes(group = compartment), fun = mean, geom = "line",
               position = position_dodge(width = 0.4)) +
  labs(x = "Days after transplanting", y = "% of reads that are microbial (not organelle)",
       colour = "Compartment",
       title = "Microbial vs rice organelle reads over time") +
  theme_bw()
print(p_org)
ggsave("results/figures/07_microbial_vs_organelle.png", p_org, width = 8, height = 5, dpi = 300)

# ---- 4. Remove non-target sequences ----
ps_clean <- prune_taxa(category == "Microbial", ps)
ps_clean <- prune_taxa(taxa_sums(ps_clean) > 0, ps_clean)
cat(sprintf("\nRemoved %d ASVs (%.1f%% of reads)\n",
            ntaxa(ps) - ntaxa(ps_clean),
            100 * (1 - sum(sample_sums(ps_clean)) / sum(sample_sums(ps)))))
ps_clean

# ---- 5. Depth and richness after cleaning ----
summ <- data.frame(meta[sample_names(ps_clean), c("sample_id", "compartment", "day")],
                   reads = sample_sums(ps_clean),
                   observed_asvs = rowSums(as(otu_table(ps_clean), "matrix") > 0))
cat("\nReads per sample after cleaning:\n"); print(summary(summ$reads))
cat("\nMean reads and observed ASVs by compartment:\n")
print(aggregate(cbind(reads, observed_asvs) ~ compartment, data = summ,
                FUN = function(x) round(mean(x))))
cat("\nLowest-depth samples:\n")
print(head(summ[order(summ$reads), ], 5), row.names = FALSE)
write.csv(summ, "results/tables/07_sample_summary.csv", row.names = FALSE)

# ---- 6. Prevalence ----
prev <- colSums(as(otu_table(ps_clean), "matrix") > 0)
cat(sprintf("\nASVs found in only one sample: %d (%.1f%%), holding %.2f%% of reads\n",
            sum(prev == 1), 100 * mean(prev == 1),
            100 * sum(taxa_sums(ps_clean)[prev == 1]) / sum(taxa_sums(ps_clean))))

saveRDS(ps_clean, "data/processed/ps_clean.rds")
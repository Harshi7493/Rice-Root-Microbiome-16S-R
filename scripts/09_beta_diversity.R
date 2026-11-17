# =============================================================
# 09. Phylogenetic tree, Faith's PD, beta diversity and PCoA
# Tree: DECIPHER alignment -> DECIPHER::DistanceMatrix (Jukes-Cantor)
#       -> FastME balanced minimum evolution with NNI (ape::fastme.bal,
#       SPR search off for speed) -> rooted on an archaeal ASV (outgroup).
# Alignment and tree are saved, so re-Sourcing skips the slow steps.
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(DECIPHER)
  library(phangorn)
  library(picante)
  library(ggplot2)
})
set.seed(100)

ps_rare <- readRDS("data/processed/ps_rare.rds")   # rarefied to 2,765 reads
md <- data.frame(sample_data(ps_rare))
comp_cols <- c("Bulk soil" = "#8C6D31", "Rhizosphere" = "#E69F00",
               "Rhizoplane" = "#009E73", "Endosphere" = "#CC79A7")
cat("ASVs in rarefied table:", ntaxa(ps_rare), "\n")

# ---- 1. Alignment (saved) ----
aln_file <- "data/processed/asv_alignment.rds"
if (file.exists(aln_file)) {
  alignment <- readRDS(aln_file)
  cat("Alignment loaded from file\n")
} else {
  t0 <- Sys.time()
  alignment <- AlignSeqs(refseq(ps_rare), anchor = NA, processors = 1, verbose = FALSE)
  saveRDS(alignment, aln_file)
  cat(sprintf("Alignment: %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ---- 2. Tree (saved) ----
tree_file <- "data/processed/asv_tree.rds"
if (file.exists(tree_file)) {
  tree <- readRDS(tree_file)
  cat("Tree loaded from file\n")
} else {
  t0 <- Sys.time()
  
  # 2a. Pairwise distances
  dm <- DistanceMatrix(alignment, correction = "Jukes-Cantor",
                       processors = 1, verbose = TRUE)
  # Very divergent pairs can give Inf/NA under Jukes-Cantor: cap them
  bad <- !is.finite(dm)
  if (any(bad)) {
    cat("Capping", sum(bad) / 2, "non-finite distances\n")
    dm[bad] <- max(dm[!bad]) * 1.5
  }
  cat(sprintf("2a. Distances done: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  # 2b. FastME tree (NNI only; SPR is far too slow for thousands of tips)
  tree <- ape::fastme.bal(as.dist(dm), nni = TRUE, spr = FALSE, tbr = FALSE)
  tree$edge.length[tree$edge.length < 0] <- 0
  rm(dm); invisible(gc())
  cat(sprintf("2b. FastME tree done: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  
  # 2c. Root on the most abundant archaeal ASV (outgroup to Bacteria)
  tax  <- as.data.frame(tax_table(ps_rare))
  arch <- rownames(tax)[tax$Kingdom %in% "Archaea"]
  if (length(arch) > 0) {
    outgroup <- arch[which.max(taxa_sums(ps_rare)[arch])]
    tree <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE)
    cat("2c. Rooted on archaeal ASV", outgroup, "\n")
  } else {
    tree <- ape::root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
    cat("2c. No archaea present; rooted on first tip\n")
  }
  tree$edge.length[is.na(tree$edge.length)] <- 0
  saveRDS(tree, tree_file)
  cat(sprintf("Tree saved: %.1f min in total\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
ps_rare <- merge_phyloseq(ps_rare, phy_tree(tree))
ape::write.tree(phy_tree(ps_rare), "results/09_asv_tree.nwk")
ps_rare

# ---- 3. Faith's phylogenetic diversity ----
otu_mat <- as(otu_table(ps_rare), "matrix")
if (taxa_are_rows(ps_rare)) otu_mat <- t(otu_mat)
faith <- pd(otu_mat, phy_tree(ps_rare), include.root = TRUE)
faith <- cbind(sample_id = rownames(faith), faith,
               md[rownames(faith), c("compartment", "day", "day_f")])

cat("\nMean Faith's PD by compartment:\n")
print(aggregate(PD ~ compartment, data = faith, FUN = function(x) round(mean(x), 1)))
cat("Kruskal-Wallis, Faith's PD ~ compartment: p =",
    signif(kruskal.test(PD ~ compartment, data = faith)$p.value, 3), "\n")
write.csv(faith, "results/tables/09_faith_pd.csv", row.names = FALSE)

# ---- 4. Distance matrices ----
cat("\nComputing beta-diversity distances...\n")
dists <- list(
  "Jaccard"            = phyloseq::distance(ps_rare, method = "jaccard", binary = TRUE),
  "Bray-Curtis"        = phyloseq::distance(ps_rare, method = "bray"),
  "Unweighted UniFrac" = phyloseq::distance(ps_rare, method = "unifrac"),
  "Weighted UniFrac"   = phyloseq::distance(ps_rare, method = "wunifrac")
)

# ---- 5. PCoA ----
pcoa_df <- do.call(rbind, lapply(names(dists), function(n) {
  o  <- ape::pcoa(dists[[n]])
  ve <- round(100 * o$values$Relative_eig[1:2], 1)
  ids <- rownames(o$vectors)
  data.frame(sample_id = ids, PC1 = o$vectors[, 1], PC2 = o$vectors[, 2],
             metric = n, PC1_pct = ve[1], PC2_pct = ve[2],
             panel = sprintf("%s (PC1 %.1f%%, PC2 %.1f%%)", n, ve[1], ve[2]),
             compartment = md[ids, "compartment"], day_f = md[ids, "day_f"])
}))
pcoa_df$panel <- factor(pcoa_df$panel, levels = unique(pcoa_df$panel))

cat("\nVariance explained by the first two PCoA axes:\n")
print(unique(pcoa_df[, c("metric", "PC1_pct", "PC2_pct")]), row.names = FALSE)

p_comp <- ggplot(pcoa_df, aes(PC1, PC2, colour = compartment, shape = day_f)) +
  geom_point(size = 2.8, alpha = 0.85) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = comp_cols) +
  labs(colour = "Compartment", shape = "Day",
       title = "Beta diversity by compartment (cf. Edwards et al. Fig. 1C)") +
  theme_bw()
print(p_comp)
ggsave("results/figures/09_pcoa_compartment.png", p_comp, width = 11, height = 8, dpi = 300)

p_day <- ggplot(pcoa_df, aes(PC1, PC2, colour = day_f, shape = compartment)) +
  geom_point(size = 2.8, alpha = 0.85) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = c("Day 1" = "#56B4E9", "Day 5" = "#0072B2",
                                 "Day 13" = "#D55E00")) +
  labs(colour = "Day", shape = "Compartment",
       title = "Beta diversity by day (cf. Edwards et al. Fig. 6C)") +
  theme_bw()
print(p_day)
ggsave("results/figures/09_pcoa_day.png", p_day, width = 11, height = 8, dpi = 300)

# ---- 6. Save for PERMANOVA ----
saveRDS(ps_rare, "data/processed/ps_rare_tree.rds")
saveRDS(dists,   "data/processed/beta_distances.rds")
cat("\nEpisode 09 complete.\n")
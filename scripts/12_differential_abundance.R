# =============================================================
# 12. Differential abundance: root compartments vs bulk soil
# ANCOM-BC2, taxa aggregated to genus (else family / order),
# model ~ compartment + day; reference = bulk soil.
# Compare with Edwards et al. (2015) Fig. 2.
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ANCOMBC)
  library(ggplot2)
})
set.seed(100)

ps  <- readRDS("data/processed/ps_clean.rds")
md  <- data.frame(sample_data(ps))
md$compartment <- factor(md$compartment,
                         levels = c("Bulk soil", "Rhizosphere", "Rhizoplane", "Endosphere"))
tax <- as.data.frame(tax_table(ps))

# ---- 1. Aggregate ASVs to genus (else family / order) ----
lab <- setNames(
  ifelse(!is.na(tax$Genus),  tax$Genus,
         ifelse(!is.na(tax$Family), paste0("f. ", tax$Family),
                ifelse(!is.na(tax$Order),  paste0("o. ", tax$Order),
                       ifelse(!is.na(tax$Class),  paste0("c. ", tax$Class), "Unclassified")))),
  rownames(tax))
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)              # taxa x samples
agg <- rowsum(otu, group = lab[rownames(otu)])     # labels x samples
phylum_of <- tapply(tax$Phylum, lab[rownames(tax)], function(x) names(which.max(table(x))))
tt <- cbind(Label = rownames(agg), Phylum = phylum_of[rownames(agg)])
rownames(tt) <- rownames(agg)

ps_g <- phyloseq(otu_table(agg, taxa_are_rows = TRUE), tax_table(tt), sample_data(md))
cat("Taxa after aggregation:", ntaxa(ps_g), "\n")

# ---- 2. ANCOM-BC2 ----
t0 <- Sys.time()
out <- ancombc2(data = ps_g, tax_level = NULL,
                fix_formula = "compartment + day_f",
                p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000,
                group = NULL, struc_zero = FALSE,
                alpha = 0.05, n_cl = 1, verbose = FALSE)
cat(sprintf("ANCOM-BC2 took %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res <- out$res
cat("Taxa tested (present in >= 10% of samples):", nrow(res), "\n")

# ---- 3. Tidy results for each compartment vs bulk soil ----
comps <- c("Rhizosphere", "Rhizoplane", "Endosphere")
tidy <- do.call(rbind, lapply(comps, function(cp) {
  s <- paste0("compartment", cp)
  data.frame(compartment = cp, taxon = res$taxon,
             phylum = tt[res$taxon, "Phylum"],
             lfc = res[[paste0("lfc_", s)]], se = res[[paste0("se_", s)]],
             q = res[[paste0("q_", s)]],
             significant = res[[paste0("diff_", s)]] & res[[paste0("passed_ss_", s)]])
}))
tidy$direction <- ifelse(!tidy$significant, "n.s.",
                         ifelse(tidy$lfc > 0, "Enriched", "Depleted"))
tidy$compartment <- factor(tidy$compartment, levels = comps)
write.csv(tidy, "results/tables/12_ancombc2_vs_bulk_soil.csv", row.names = FALSE)

# ---- 4. Counts of enriched / depleted taxa (cf. Fig. 2B-C) ----
counts <- table(tidy$compartment, tidy$direction)[, c("Enriched", "Depleted")]
cat("\nSignificant taxa vs bulk soil (Holm q < 0.05, robust to pseudo-counts):\n")
print(counts)

enr <- lapply(comps, function(cp) tidy$taxon[tidy$compartment == cp & tidy$direction == "Enriched"])
dep <- lapply(comps, function(cp) tidy$taxon[tidy$compartment == cp & tidy$direction == "Depleted"])
names(enr) <- names(dep) <- comps
cat("\nEnriched in all three root compartments:", length(Reduce(intersect, enr)), "\n")
cat("Depleted in all three root compartments:", length(Reduce(intersect, dep)), "\n")
cat("Rhizoplane & endosphere share", length(intersect(enr$Rhizoplane, enr$Endosphere)),
    "enriched taxa\n")

core_enr <- tidy[tidy$compartment == "Endosphere" & tidy$taxon %in% Reduce(intersect, enr),
                 c("taxon", "phylum", "lfc", "q")]
cat("\nTaxa enriched in ALL root compartments (endosphere LFC shown):\n")
print(core_enr[order(-core_enr$lfc), ], row.names = FALSE)

cat("\nTop 10 enriched in the endosphere:\n")
e <- tidy[tidy$compartment == "Endosphere" & tidy$direction == "Enriched", ]
print(head(e[order(-e$lfc), c("taxon", "phylum", "lfc", "q")], 10), row.names = FALSE)
cat("\nTop 10 depleted in the endosphere:\n")
d <- tidy[tidy$compartment == "Endosphere" & tidy$direction == "Depleted", ]
print(head(d[order(d$lfc), c("taxon", "phylum", "lfc", "q")], 10), row.names = FALSE)

cat("\nSignificant taxa by phylum, endosphere:\n")
print(table(droplevels(factor(e$phylum))), max = 20)
print(table(droplevels(factor(d$phylum))), max = 20)

# ---- 5. Plots ----
cnt_df <- as.data.frame(counts); names(cnt_df) <- c("compartment", "direction", "n")
p_counts <- ggplot(cnt_df, aes(compartment, ifelse(direction == "Depleted", -n, n),
                               fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = n, vjust = ifelse(direction == "Depleted", 1.4, -0.4))) +
  scale_fill_manual(values = c(Enriched = "#0072B2", Depleted = "#D55E00")) +
  labs(x = NULL, y = "Number of taxa (depleted shown below zero)", fill = NULL,
       title = "Taxa enriched / depleted vs bulk soil (cf. Edwards et al. Fig. 2)") +
  theme_bw()
print(p_counts)
ggsave("results/figures/12_counts_vs_soil.png", p_counts, width = 7, height = 5, dpi = 300)

# Heatmap: 40 taxa with the largest significant |LFC| in any compartment
sig <- tidy[tidy$significant, ]
top <- head(names(sort(tapply(abs(sig$lfc), sig$taxon, max), decreasing = TRUE)), 40)
hm <- tidy[tidy$taxon %in% top, ]
hm$lfc_plot <- ifelse(hm$significant, hm$lfc, NA)
ord <- with(hm[hm$compartment == "Endosphere", ], taxon[order(lfc)])
hm$taxon <- factor(hm$taxon, levels = ord)
p_hm <- ggplot(hm, aes(compartment, taxon, fill = lfc_plot)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "#D55E00", mid = "white", high = "#0072B2",
                       midpoint = 0, na.value = "grey92", name = "LFC vs\nbulk soil") +
  labs(x = NULL, y = NULL, title = "Largest changes vs bulk soil (grey = not significant)") +
  theme_minimal() + theme(axis.text.y = element_text(size = 7, face = "italic"))
print(p_hm)
ggsave("results/figures/12_heatmap_vs_soil.png", p_hm, width = 7, height = 9, dpi = 300)
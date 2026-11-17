# =============================================================
# 11. Taxonomic composition by compartment and day
# Relative abundance from the unrarefied, organelle-free table.
# Compare with Edwards et al. (2015) Fig. 1D and Results.
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
  library(RColorBrewer)
})

ps  <- readRDS("data/processed/ps_clean.rds")
md  <- data.frame(sample_data(ps))
tax <- as.data.frame(tax_table(ps))
rel <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) rel <- t(rel)
rel <- rel / rowSums(rel)                       # samples x ASVs, proportions

# ---- Label lookups (named by ASV) ----
genus_lab <- setNames(
  ifelse(!is.na(tax$Genus),  tax$Genus,
         ifelse(!is.na(tax$Family), paste0("f. ", tax$Family),
                ifelse(!is.na(tax$Order),  paste0("o. ", tax$Order), "Unclassified"))),
  rownames(tax))
phylum_lab <- setNames(ifelse(is.na(tax$Phylum), "Unclassified", tax$Phylum),
                       rownames(tax))
stopifnot(all(colnames(rel) %in% names(genus_lab)))   # every ASV has a label

# Sample order: compartment, then day
md <- md[order(md$compartment, md$day), ]
sample_order <- md$sample_id

# ---- Helper: sum ASVs to a rank, keep top N (by max compartment mean) ----
collapse_rank <- function(labels, top_n) {
  m <- t(rowsum(t(rel), group = labels[colnames(rel)]))   # samples x taxa
  m <- m[sample_order, , drop = FALSE]
  comp_means <- rowsum(m, md[rownames(m), "compartment"]) /
    as.vector(table(md[rownames(m), "compartment"]))
  top <- names(sort(apply(comp_means, 2, max), decreasing = TRUE))
  top <- head(setdiff(top, "Unclassified"), top_n)
  other <- rowSums(m[, !colnames(m) %in% top, drop = FALSE])
  m <- cbind(m[, top, drop = FALSE], Other = other)
  data.frame(sample_id   = rep(rownames(m), times = ncol(m)),
             Taxon       = factor(rep(colnames(m), each = nrow(m)), levels = colnames(m)),
             Abundance   = as.vector(m),
             compartment = rep(md[rownames(m), "compartment"], times = ncol(m)),
             day_f       = rep(md[rownames(m), "day_f"], times = ncol(m)))
}
phylum_df <- collapse_rank(phylum_lab, 10)
genus_df  <- collapse_rank(genus_lab, 15)

make_pal <- function(levels) {
  base <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
  c(setNames(base[seq_len(length(levels) - 1)], head(levels, -1)), Other = "grey80")
}

# ---- 1. Per-sample bar plots ----
bar_plot <- function(d, title) {
  d$sample_id <- factor(d$sample_id, levels = sample_order)
  ggplot(d, aes(sample_id, Abundance, fill = Taxon)) +
    geom_col(width = 0.95) +
    facet_grid(~ compartment, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = make_pal(levels(d$Taxon))) +
    scale_x_discrete(labels = function(x) as.character(md[x, "day_f"])) +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = "Samples (labelled by day)", y = "Relative abundance", fill = NULL, title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          legend.text = element_text(size = 8))
}
p_phylum <- bar_plot(phylum_df, "Phylum composition per sample")
p_genus  <- bar_plot(genus_df, "Genus composition per sample (top 15)") +
  theme(legend.text = element_text(size = 8, face = "italic"))
print(p_phylum); print(p_genus)
ggsave("results/figures/11_phylum_per_sample.png", p_phylum, width = 13, height = 6, dpi = 300)
ggsave("results/figures/11_genus_per_sample.png",  p_genus,  width = 13, height = 7, dpi = 300)

# ---- 2. Mean composition by compartment x day ----
mean_by <- function(d) aggregate(Abundance ~ compartment + day_f + Taxon, data = d, FUN = mean)
p_phylum_time <- ggplot(mean_by(phylum_df), aes(day_f, Abundance, fill = Taxon)) +
  geom_col(width = 0.8) +
  facet_grid(~ compartment) +
  scale_fill_manual(values = make_pal(levels(phylum_df$Taxon))) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = NULL, y = "Mean relative abundance", fill = NULL,
       title = "Mean phylum composition by compartment and day (cf. Edwards et al. Fig. 1D)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_phylum_time)
ggsave("results/figures/11_phylum_by_compartment_day.png", p_phylum_time,
       width = 11, height = 5.5, dpi = 300)

# ---- 3. Summary tables (mean % by compartment) ----
to_wide <- function(d) {
  w <- round(100 * tapply(d$Abundance, list(d$Taxon, d$compartment), mean), 1)
  w[order(-rowMeans(w)), ]
}
phylum_tab <- to_wide(phylum_df); genus_tab <- to_wide(genus_df)
cat("Mean phylum abundance (%) by compartment:\n"); print(phylum_tab)
cat("\nMean genus abundance (%) by compartment (top 15 + Other):\n"); print(genus_tab)
write.csv(phylum_tab, "results/tables/11_phylum_by_compartment.csv")
write.csv(genus_tab,  "results/tables/11_genus_by_compartment.csv")

# ---- 4. Taxa highlighted by Edwards et al. (2015) ----
group_pct <- function(keep) {
  s <- rowSums(rel[, keep, drop = FALSE])
  round(100 * tapply(s, md[names(s), "compartment"], mean), 2)
}
highlight <- rbind(
  "Acidobacteriota (phylum)"      = group_pct(tax$Phylum %in% "Acidobacteriota"),
  "Planctomycetota (phylum)"      = group_pct(tax$Phylum %in% "Planctomycetota"),
  "Gemmatimonadota (phylum)"      = group_pct(tax$Phylum %in% "Gemmatimonadota"),
  "Pseudomonadota (phylum)"       = group_pct(tax$Phylum %in% "Pseudomonadota"),
  "Spirochaetota (phylum)"        = group_pct(tax$Phylum %in% "Spirochaetota"),
  "Fibrobacterota (phylum)"       = group_pct(tax$Phylum %in% "Fibrobacterota"),
  "Rhodocyclaceae (family)"       = group_pct(tax$Family %in% "Rhodocyclaceae"),
  "Comamonadaceae (family)"       = group_pct(tax$Family %in% "Comamonadaceae"),
  "Pleomorphomonas (genus)"       = group_pct(tax$Genus  %in% "Pleomorphomonas"),
  "Anaeromyxobacter (genus)"      = group_pct(tax$Genus  %in% "Anaeromyxobacter"),
  "Methanogens (archaea, genera)" = group_pct(grepl("^Methano", tax$Genus))
)
cat("\nTaxa highlighted by Edwards et al. (2015): mean % by compartment\n")
print(highlight)
write.csv(highlight, "results/tables/11_edwards_highlight_taxa.csv")

# ---- 5. Endosphere over time (top phyla) ----
cat("\nMean phylum % over time, endosphere:\n")
e <- mean_by(phylum_df[phylum_df$compartment == "Endosphere", ])
print(round(100 * tapply(e$Abundance, list(e$Taxon, e$day_f), sum), 1))
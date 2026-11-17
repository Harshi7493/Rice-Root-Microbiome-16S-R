# =============================================================
# 13. Publication-quality multi-panel figure (Figure 1)
# Reads saved results from Episodes 07-12; only panel D distances
# are recomputed (quick). Legends shown once (from panels C and E).
# =============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(phyloseq)
})

comp_levels <- c("Bulk soil", "Rhizosphere", "Rhizoplane", "Endosphere")
day_levels  <- c("Day 1", "Day 5", "Day 13")
comp_cols   <- c("Bulk soil" = "#8C6D31", "Rhizosphere" = "#E69F00",
                 "Rhizoplane" = "#009E73", "Endosphere" = "#CC79A7")
day_shapes  <- c("Day 1" = 16, "Day 5" = 17, "Day 13" = 15)
fix <- function(d) {
  d$compartment <- factor(d$compartment, levels = comp_levels)
  if ("day_f" %in% names(d)) d$day_f <- factor(d$day_f, levels = day_levels)
  d
}
theme_pub <- theme_classic(base_size = 9) +
  theme(plot.title    = element_text(face = "bold", size = 10),
        plot.subtitle = element_text(size = 7.5, lineheight = 1.1),
        legend.title  = element_text(size = 8, face = "bold"),
        legend.text   = element_text(size = 8))
dodge <- position_dodge(width = 0.35)

# ---- A. Colonisation of the root interior (no legend) ----
org <- fix(read.csv("results/tables/07_organelle_by_sample.csv", check.names = FALSE))
endo <- org[org$compartment == "Endosphere", ]
em <- tapply(endo$Microbial, endo$day_f, mean)
pA <- ggplot(org, aes(day_f, Microbial, colour = compartment, group = compartment)) +
  geom_point(size = 1.6, alpha = 0.7, position = dodge) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = dodge) +
  scale_colour_manual(values = comp_cols, guide = "none") +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Days after transplanting", y = "Microbial reads (% of total)",
       title = "Root colonisation",
       subtitle = sprintf("Endosphere: %.0f%% \u2192 %.0f%% \u2192 %.0f%% microbial reads\n(remainder = rice organelle DNA)",
                          em[1], em[2], em[3])) +
  theme_pub

# ---- B. Shannon diversity (no legend) ----
alpha <- fix(read.csv("results/tables/08_alpha_diversity.csv"))
an <- read.csv("results/tables/08_alpha_anova.csv")
p_comp <- an$p[an$metric == "Shannon" & an$term == "compartment"]
pB <- ggplot(alpha, aes(compartment, Shannon)) +
  geom_boxplot(aes(fill = compartment), alpha = 0.4, width = 0.6,
               outlier.shape = NA, linewidth = 0.3) +
  geom_point(aes(colour = compartment, shape = day_f),
             position = position_jitter(width = 0.12, seed = 1), size = 1.6) +
  scale_fill_manual(values = comp_cols, guide = "none") +
  scale_colour_manual(values = comp_cols, guide = "none") +
  scale_shape_manual(values = day_shapes, guide = "none") +
  labs(x = NULL, y = "Shannon diversity", title = "Alpha diversity",
       subtitle = sprintf("Rarefied to 2,765 reads\nANOVA (compartment) p = %.1e", p_comp)) +
  theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ---- C. PCoA, weighted UniFrac (carries the Compartment and Day legends) ----
ps    <- readRDS("data/processed/ps_rare_tree.rds")
md    <- fix(data.frame(sample_data(ps)))
dists <- readRDS("data/processed/beta_distances.rds")
pc    <- ape::pcoa(dists[["Weighted UniFrac"]])
ve    <- 100 * pc$values$Relative_eig[1:2]
pc_df <- data.frame(PC1 = pc$vectors[, 1], PC2 = pc$vectors[, 2],
                    md[rownames(pc$vectors), c("compartment", "day_f")])
perm <- read.csv("results/tables/10_permanova_main.csv")
r2 <- perm$R2[perm$metric == "Weighted UniFrac" & perm$term == "compartment"]
pC <- ggplot(pc_df, aes(PC1, PC2, colour = compartment, shape = day_f)) +
  geom_point(size = 2, alpha = 0.85) +
  scale_colour_manual(values = comp_cols, name = "Compartment", drop = FALSE,
                      guide = guide_legend(order = 1, nrow = 1,
                                           override.aes = list(size = 3))) +
  scale_shape_manual(values = day_shapes, name = "Day",
                     guide = guide_legend(order = 2, nrow = 1,
                                          override.aes = list(size = 3))) +
  labs(x = sprintf("PCoA1 (%.1f%%)", ve[1]), y = sprintf("PCoA2 (%.1f%%)", ve[2]),
       title = "Beta diversity (weighted UniFrac)",
       subtitle = sprintf("PERMANOVA (compartment)\nR\u00b2 = %.2f, p = 0.001", r2)) +
  theme_pub

# ---- D. Divergence from same-day bulk soil (no legend) ----
bc <- as.matrix(dists[["Bray-Curtis"]])
div <- do.call(rbind, lapply(rownames(md)[md$compartment != "Bulk soil"], function(s) {
  soil <- rownames(md)[md$compartment == "Bulk soil" & md$day == md[s, "day"]]
  data.frame(compartment = md[s, "compartment"], day_f = md[s, "day_f"],
             dist = mean(bc[s, soil]))
}))
pD <- ggplot(div, aes(day_f, dist, colour = compartment, group = compartment)) +
  geom_point(size = 1.6, alpha = 0.7, position = dodge) +
  stat_summary(fun = mean, geom = "line", linewidth = 0.8, position = dodge) +
  scale_colour_manual(values = comp_cols, guide = "none") +
  labs(x = "Days after transplanting", y = "Bray-Curtis distance to soil",
       title = "Divergence from soil",
       subtitle = "Distance to same-day bulk soil\nSpearman \u03c1 = 0.74\u20130.79, BH p \u2264 0.031") +
  theme_pub

# ---- E. Enriched / depleted taxa vs bulk soil (carries its own legend) ----
da <- read.csv("results/tables/12_ancombc2_vs_bulk_soil.csv")
da <- da[da$direction != "n.s.", ]
cnt <- as.data.frame(table(factor(da$compartment, levels = comp_levels[-1]),
                           factor(da$direction, levels = c("Enriched", "Depleted"))))
names(cnt) <- c("compartment", "direction", "n")
cnt$y <- ifelse(cnt$direction == "Depleted", -cnt$n, cnt$n)
pE <- ggplot(cnt, aes(compartment, y, fill = direction)) +
  geom_col(width = 0.55) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(aes(label = n, vjust = ifelse(direction == "Depleted", 1.3, -0.3)), size = 2.8) +
  scale_fill_manual(values = c(Enriched = "#0072B2", Depleted = "#D55E00"),
                    name = "Vs bulk soil", guide = guide_legend(order = 3, nrow = 1)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(x = NULL, y = "Number of taxa", title = "Differential abundance",
       subtitle = "ANCOM-BC2 vs bulk soil\ngenus level, Holm q < 0.05") +
  theme_pub

# ---- Assemble: 3 panels on top, 2 wider panels below, one legend row ----
design <- "AABBCC
           DDDEEE"
fig <- pA + pB + pC + pD + pE +
  plot_layout(design = design, guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Assembly of the rice root microbiome over 13 days (reanalysis of Edwards et al. 2015)",
                  theme = theme(plot.title = element_text(face = "bold", size = 12))) &
  theme(plot.tag = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.box = "horizontal")

print(fig)
ggsave("results/figures/13_figure1.png", fig, width = 270, height = 200, units = "mm",
       dpi = 300, bg = "white")
ggsave("results/figures/13_figure1.pdf", fig, width = 270, height = 200, units = "mm")
cat("Saved results/figures/13_figure1.png and .pdf\n")
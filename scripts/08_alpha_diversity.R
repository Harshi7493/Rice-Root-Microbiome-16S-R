# =============================================================
# 08. Alpha diversity: compartment and time
# Main analysis rarefied to 2,765 reads (all 36 samples);
# sensitivity analysis at 10,000 reads (samples with enough depth).
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})

ps_clean   <- readRDS("data/processed/ps_clean.rds")
depth_main <- min(sample_sums(ps_clean))   # 2,765 (lowest sample)
depth_sens <- 10000
comp_cols  <- c("Bulk soil" = "#8C6D31", "Rhizosphere" = "#E69F00",
                "Rhizoplane" = "#009E73", "Endosphere" = "#CC79A7")

# ---- 1. Rarefaction curves ----
otu <- as(otu_table(ps_clean), "matrix")
rc  <- rarecurve(otu, step = 1000, tidy = TRUE)
names(rc) <- c("sample_id", "reads", "ASVs")
md <- data.frame(sample_data(ps_clean))
rc$compartment <- md[rc$sample_id, "compartment"]

p_rare <- ggplot(rc, aes(reads, ASVs, group = sample_id, colour = compartment)) +
  geom_line(alpha = 0.8) +
  geom_vline(xintercept = c(depth_main, depth_sens), linetype = c("dashed", "dotted")) +
  scale_colour_manual(values = comp_cols) +
  labs(x = "Reads sampled", y = "Observed ASVs", colour = "Compartment",
       title = "Rarefaction curves (dashed = 2,765; dotted = 10,000)") +
  theme_bw()
print(p_rare)
ggsave("results/figures/08_rarefaction_curves.png", p_rare, width = 9, height = 5, dpi = 300)

# ---- 2. Alpha metrics on a rarefied table ----
alpha_at <- function(depth) {
  ps_r <- rarefy_even_depth(ps_clean, sample.size = depth, rngseed = 100,
                            replace = FALSE, trimOTUs = TRUE, verbose = FALSE)
  a <- estimate_richness(ps_r, measures = c("Observed", "Shannon"))
  a$Pielou <- a$Shannon / log(a$Observed)
  a$sample_id <- rownames(a)
  m <- data.frame(sample_data(ps_r))[a$sample_id, c("compartment", "day", "day_f")]
  list(ps = ps_r, alpha = cbind(a, m))
}

main  <- alpha_at(depth_main)
alpha <- main$alpha
cat("Main analysis: rarefied to", depth_main, "reads,", nrow(alpha), "samples\n")

cat("\nMean alpha diversity by compartment:\n")
print(aggregate(cbind(Observed, Shannon, Pielou) ~ compartment, data = alpha,
                FUN = function(x) round(mean(x), 2)))
cat("\nMean Shannon by compartment and day:\n")
print(round(tapply(alpha$Shannon, list(alpha$compartment, alpha$day_f), mean), 2))
write.csv(alpha, "results/tables/08_alpha_diversity.csv", row.names = FALSE)

# ---- 3. Statistics ----
metrics <- c("Observed", "Shannon", "Pielou")

# (a) Two-way ANOVA: compartment, day and their interaction
anova_tab <- do.call(rbind, lapply(metrics, function(m) {
  a <- summary(aov(alpha[[m]] ~ compartment * day_f, data = alpha))[[1]]
  data.frame(metric = m, term = trimws(rownames(a))[1:3],
             F = round(a$`F value`[1:3], 2), p = signif(a$`Pr(>F)`[1:3], 3))
}))
cat("\nTwo-way ANOVA (compartment x day):\n"); print(anova_tab, row.names = FALSE)

# (b) Kruskal-Wallis by compartment (rank-based check)
kw <- sapply(metrics, function(m) kruskal.test(alpha[[m]] ~ alpha$compartment)$p.value)
cat("\nKruskal-Wallis by compartment, p-values:\n"); print(signif(kw, 3))

# (c) Pairwise compartments (Shannon, BH-adjusted)
cat("\nPairwise Wilcoxon, Shannon by compartment (BH-adjusted):\n")
print(pairwise.wilcox.test(alpha$Shannon, alpha$compartment,
                           p.adjust.method = "BH", exact = FALSE)$p.value, digits = 3)

# (d) Effect of day within each compartment (Shannon, one-way ANOVA)
day_within <- do.call(rbind, lapply(levels(alpha$compartment), function(cp) {
  d <- alpha[alpha$compartment == cp, ]
  a <- summary(aov(Shannon ~ day_f, data = d))[[1]]
  data.frame(compartment = cp, F = round(a$`F value`[1], 2), p = signif(a$`Pr(>F)`[1], 3))
}))
cat("\nEffect of day on Shannon within each compartment:\n")
print(day_within, row.names = FALSE)

write.csv(anova_tab,  "results/tables/08_alpha_anova.csv", row.names = FALSE)
write.csv(day_within, "results/tables/08_alpha_day_within.csv", row.names = FALSE)

# ---- 4. Plots ----
alpha_long <- do.call(rbind, lapply(metrics, function(m)
  data.frame(alpha[, c("sample_id", "compartment", "day_f")], metric = m, value = alpha[[m]])))
alpha_long$metric <- factor(alpha_long$metric, levels = metrics,
                            labels = c("Observed ASVs", "Shannon", "Pielou evenness"))

p_alpha <- ggplot(alpha_long, aes(compartment, value)) +
  geom_boxplot(aes(fill = compartment), alpha = 0.4, outlier.shape = NA) +
  geom_point(aes(shape = day_f), position = position_jitter(width = 0.15, seed = 1),
             size = 2.2) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = comp_cols, guide = "none") +
  labs(x = NULL, y = NULL, shape = "Day",
       title = paste0("Alpha diversity by compartment (rarefied to ", depth_main, " reads)")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
print(p_alpha)
ggsave("results/figures/08_alpha_by_compartment.png", p_alpha, width = 11, height = 4.5, dpi = 300)

p_time <- ggplot(alpha, aes(day_f, Shannon, colour = compartment, group = compartment)) +
  geom_point(size = 2.5, alpha = 0.7, position = position_dodge(width = 0.3)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1,
               position = position_dodge(width = 0.3)) +
  scale_colour_manual(values = comp_cols) +
  labs(x = "Days after transplanting", y = "Shannon diversity", colour = "Compartment",
       title = "Shannon diversity over time") +
  theme_bw()
print(p_time)
ggsave("results/figures/08_shannon_over_time.png", p_time, width = 8, height = 5, dpi = 300)

# ---- 5. Sensitivity analysis at 10,000 reads ----
sens <- alpha_at(depth_sens)
cat("\nSensitivity at", depth_sens, "reads:", nrow(sens$alpha), "samples kept; dropped:",
    setdiff(alpha$sample_id, sens$alpha$sample_id), "\n")
print(aggregate(cbind(Observed, Shannon) ~ compartment, data = sens$alpha,
                FUN = function(x) round(mean(x), 2)))
cat("Kruskal-Wallis by compartment, p-values:\n")
print(signif(sapply(metrics, function(m)
  kruskal.test(sens$alpha[[m]] ~ sens$alpha$compartment)$p.value), 3))

# ---- 6. Save rarefied object for beta diversity ----
saveRDS(main$ps, "data/processed/ps_rare.rds")
# =============================================================
# 10. PERMANOVA: compartment, day, divergence from soil, dispersion
# =============================================================

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
})
set.seed(100)

ps    <- readRDS("data/processed/ps_rare_tree.rds")
dists <- readRDS("data/processed/beta_distances.rds")
md    <- data.frame(sample_data(ps))
comp_cols <- c("Bulk soil" = "#8C6D31", "Rhizosphere" = "#E69F00",
               "Rhizoplane" = "#009E73", "Endosphere" = "#CC79A7")
sub_dist <- function(d, keep) as.dist(as.matrix(d)[keep, keep])

# ---- 1. Two-factor PERMANOVA (balanced design: sequential terms) ----
perm_main <- do.call(rbind, lapply(names(dists), function(n) {
  d <- dists[[n]]; m <- md[labels(d), ]
  r <- adonis2(d ~ compartment * day_f, data = m, permutations = 999, by = "terms")
  data.frame(metric = n, term = rownames(r)[1:3], R2 = round(r$R2[1:3], 3),
             F = round(r$F[1:3], 2), p = r$`Pr(>F)`[1:3])
}))
cat("PERMANOVA: compartment x day\n"); print(perm_main, row.names = FALSE)
write.csv(perm_main, "results/tables/10_permanova_main.csv", row.names = FALSE)

# ---- 2. Pairwise compartments (BH-adjusted) ----
pairwise_perm <- function(d) {
  m <- md[labels(d), ]
  prs <- combn(levels(m$compartment), 2)
  out <- do.call(rbind, lapply(seq_len(ncol(prs)), function(i) {
    keep <- m$compartment %in% prs[, i]
    r <- adonis2(sub_dist(d, keep) ~ compartment, data = droplevels(m[keep, ]),
                 permutations = 999)
    data.frame(group1 = prs[1, i], group2 = prs[2, i],
               R2 = round(r$R2[1], 3), p = r$`Pr(>F)`[1])
  }))
  out$p_adj <- round(p.adjust(out$p, method = "BH"), 4)
  out
}
pw <- do.call(rbind, lapply(c("Bray-Curtis", "Weighted UniFrac"), function(n)
  cbind(metric = n, pairwise_perm(dists[[n]]))))
cat("\nPairwise PERMANOVA between compartments:\n"); print(pw, row.names = FALSE)
write.csv(pw, "results/tables/10_pairwise_compartments.csv", row.names = FALSE)

# ---- 3. Effect of day within each compartment ----
day_within <- do.call(rbind, lapply(names(dists), function(n) {
  do.call(rbind, lapply(levels(md$compartment), function(cp) {
    keep <- md[labels(dists[[n]]), "compartment"] == cp
    m <- md[labels(dists[[n]]), ][keep, ]
    r <- adonis2(sub_dist(dists[[n]], keep) ~ day_f, data = m, permutations = 999)
    data.frame(metric = n, compartment = cp, R2 = round(r$R2[1], 3), p = r$`Pr(>F)`[1])
  }))
}))
day_within$p_adj <- round(p.adjust(day_within$p, method = "BH"), 4)
cat("\nPERMANOVA of day within each compartment (BH across all tests):\n")
print(day_within, row.names = FALSE)
write.csv(day_within, "results/tables/10_day_within_compartment.csv", row.names = FALSE)

# ---- 4. Divergence from same-day bulk soil ----
divergence <- do.call(rbind, lapply(names(dists), function(n) {
  dm <- as.matrix(dists[[n]]); m <- md[rownames(dm), ]
  root <- rownames(m)[m$compartment != "Bulk soil"]
  do.call(rbind, lapply(root, function(s) {
    soil <- rownames(m)[m$compartment == "Bulk soil" & m$day == m[s, "day"]]
    data.frame(metric = n, sample_id = s, compartment = m[s, "compartment"],
               day = m[s, "day"], day_f = m[s, "day_f"],
               dist_to_soil = mean(dm[s, soil]))
  }))
}))
divergence$compartment <- droplevels(divergence$compartment)

div_tests <- do.call(rbind, lapply(split(divergence,
                                         list(divergence$metric, divergence$compartment)), function(d) {
                                           ct <- cor.test(d$dist_to_soil, d$day, method = "spearman", exact = FALSE)
                                           data.frame(metric = d$metric[1], compartment = as.character(d$compartment[1]),
                                                      day1 = round(mean(d$dist_to_soil[d$day == 1]), 3),
                                                      day5 = round(mean(d$dist_to_soil[d$day == 5]), 3),
                                                      day13 = round(mean(d$dist_to_soil[d$day == 13]), 3),
                                                      spearman_rho = round(unname(ct$estimate), 2), p = signif(ct$p.value, 3))
                                         }))
div_tests$p_adj <- signif(p.adjust(div_tests$p, method = "BH"), 3)
div_tests <- div_tests[order(div_tests$metric, div_tests$compartment), ]
cat("\nMean distance to same-day bulk soil, and trend over time:\n")
print(div_tests, row.names = FALSE)
write.csv(div_tests, "results/tables/10_divergence_from_soil.csv", row.names = FALSE)

p_div <- ggplot(divergence, aes(day_f, dist_to_soil, colour = compartment, group = compartment)) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.3)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1,
               position = position_dodge(width = 0.3)) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_colour_manual(values = comp_cols) +
  labs(x = "Days after transplanting", y = "Mean distance to same-day bulk soil",
       colour = "Compartment", title = "Divergence of root communities from soil over time") +
  theme_bw()
print(p_div)
ggsave("results/figures/10_divergence_from_soil.png", p_div, width = 10, height = 7, dpi = 300)

# ---- 5. Dispersion by compartment ----
disp <- do.call(rbind, lapply(names(dists), function(n) {
  d  <- dists[[n]]
  bd <- betadisper(d, md[labels(d), "compartment"])
  pt <- permutest(bd, permutations = 999)
  cbind(data.frame(metric = n, F = round(pt$tab$F[1], 2), p = pt$tab$`Pr(>F)`[1]),
        t(round(tapply(bd$distances, bd$group, mean), 3)))
}))
cat("\nDispersion (betadisper) by compartment, mean distance to centroid:\n")
print(disp, row.names = FALSE)
write.csv(disp, "results/tables/10_betadisper.csv", row.names = FALSE)
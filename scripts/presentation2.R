# Run generative_simulation.R in bcmslab/rethink
# Run all manuscript.Rmd in bcmslab/task

# Summary of the observations

# presentation2/figures/plot_expression_response_correlations2.png ----
png(filename = 'presentation2/figures/plot_expression_response_correlations2.png',
    height = 4, width = 7, units = 'in', res = 300)
par(mfrow = c(1, 2))
plot(density(expression_response_correlation$cor),
     main = 'Expression-response\n Correlation', 
     xlab = 'Spearman Rank Coefficient',
     lwd = 2)
abline(v = 0, lty = 2, col = 'red')

with(correlations_by_quantiles,
     boxplot(score_ave ~ quant_label,
             lwd = 2, xaxt = 'n',
             main = 'Response to pertubation\n by Quantiles',
             xlab = 'Expression (Quanitles)',
             ylab = 'Response (fold-change)'))
axis(1, at=1:5, labels= (1:5) * 20)

abline(h = 0, lty = 2, col = 'red')
dev.off()

# presentation2/figures/plot_archetypes2.png ---- 
rss <- rowMeans(rss(cell_archetypes), na.rm = TRUE)

simp <- simplexplot(cell_archtype_6)

png(filename = 'presentation2/figures/plot_archetypes2.png',
    height = 4, width = 7, units = 'in', res = 300)

par(mfrow = c(1,2))
plot(rss,
     type = 'l', lwd = 2,
     xlab = 'Number of Archetypes (K)',
     ylab = 'RSS',
     main = 'Residual Sum of Squares\n at each K')
abline(h = rss[6], lty = 2, col = 'red')
abline(v = 6, lty = 2, col = 'red')


plot(simp$proj_h[, 1],
     simp$proj_h[, 2],
     lwd = 2,
     xlim = c(-10, 10),
     ylim = c(-10, 10),
     xlab = 'Projection (x)',
     ylab = 'Projection (y)',
     main = 'Archetypes Simplex (2D)\n (K = 6)')
text(simp$proj_labels[,1],
     simp$proj_labels[, 2],
     labels = paste0('A', 1:6), 
     col = 'red')
dev.off()

# presentation2/figures/archetype_variance.png ----
png(filename = 'presentation2/figures/archetype_variance.png',
    height = 4, width = 7, units = 'in', res = 300)
par(mfrow = c(1,2))
plot(mds$D1,
     mds$D2,
     lwd = 2,
     col = as.factor(mds$archetype),
     xlab = 'D1', ylab = 'D2',
     main = 'Multi-dimensional Scaling (MDS)\n Color by nearest archetype')
var_part_df <- as_tibble(melt(var_part))
var_part_df$variable <- str_to_title(var_part_df$variable)

boxplot(var_part_df$value ~ var_part_df$variable,
        lwd = 2, xaxt = 'n',
        xlab = 'Variable', 
        ylab = 'Intra-class Correlation (ICC)',
        main = 'Fraction of Variance Explained\n by variables')

axis(1,cex.axis=.75, at = 1:3, labels = c('Archetype', 'Residuals', 'Tissue'))

dev.off()

# presentation2/figures/archetype_distance_enrichment.png ----
enrich_size <- with(archetype_enrichment, split(ifelse(pvalue < .05, NES, 0), ID))
enrich_shape <- lapply(with(archetype_enrichment, split(NES > 0, ID)), as.numeric)
enrich_shape <- with(archetype_enrichment, split(ifelse(NES > 0, 2, 6), ID))
xn <- length(unique(archetype_enrichment$archetype))
yn <- length(unique(archetype_enrichment$ID))

png(filename = 'presentation2/figures/archetype_distance_enrichment.png',
    height = 4, width = 7, units = 'in', res = 300)
par(mfrow = c(1,2))

plot(NULL,
     type = 'n',
     xlim = c(1,xn),
     ylim = c(1,yn),
     yaxt="n", xaxt="n",
     xlab = 'Archetype',
     ylab = 'Cancer Hallmark',
     main = 'Enrichment in Archetypes\n Enrichment Score (NES)')

for (i in seq_len(yn)) {
  points(1:xn, rep(i, xn),
         cex = abs(enrich_size[[i]]),
         pch = enrich_shape[[i]],
         lwd = 2)
}
axis(1, at = 1:xn, labels = paste0('A', 1:xn))
axis(2, at = 1:yn, labels = paste0('#', 1:yn), las = 2)

cell_distances_mat <- acast(cell_distances, cell_id ~ archetype, value.var = 'dss')
mat2 <- t(cell_distances_mat)
hc <- hclust(dist(t(mat2)))

mat2 <- mat2[, hc$order]
rownames(mat2) <- paste0('A', 1:nrow(mat2))


brca_cells <- unique(expression_response_growth$cell_id)
par(mar = c(5, 5, 4, 1))
image(mat2, 
      yaxt="n", xaxt="n",
      col = gray.colors(3),
      main = 'Distance Between Cells\n Sum of Squar. Diff. (SSD)',
      xlab = 'Archetype',
      ylab =  paste('Cells (Breast Cancer)', paste(brca_cells, collapse = '/'), sep = '\n'))

axis(1, at = seq(0, 1, .2), paste0('A', 1:nrow(mat2)))
axis(2, at = which(colnames(mat2) %in% brca_cells)/ncol(mat2),
     labels = rep('',4), cex.axis = .1)

dev.off()

# presentation2/figures/archetype_correlation_effect.png ----
# show enrichment in archetype correlate with enrichment in response to drugs
arch_enrich_groups <- archetype_enrichment_distance %>%
  filter(dss < 120) %>%
  group_by(cell_id, archetype, ID) %>%
  summarise_at(vars(ends_with('NES')), mean) %>%
  ungroup() %>%
  mutate(group = cut(abs(R_NES),
                     breaks = quantile(abs(R_NES), probs = seq(0, 1, by = 0.2)),
                     labels = c("20","40","60","80","100"),
                     include.lowest=TRUE)) %>%
  group_by(group) %>%
  summarise(corr = cor(A_NES, R_NES, use = 'complete'))

png(filename = 'presentation2/figures/archetype_correlation_effect.png',
    height = 4, width = 7, units = 'in', res = 300)
par(mar = c(5, 5, 4, 2), mfrow = c(1, 2))

with(arch_enrich_groups,
     barplot(corr ~ group,
             ylim = c(-.35, 0),
             yaxt = 'n',
             main = 'Drug Effect on Specialized\n Hallmark Enrichment',
             ylab = 'Correlation of Enrichment Scores\n (archetype and drug response)',
             xlab = 'Enrichment in Nearest Archetype\n (Absolute value; Quantile)'))
box()
axis(1, at = c(.75, 1.8, 3, 4.3, 5.5), lwd.ticks = 1, labels = NA)
axis(2, at = c(0, -.1, -.2, -.3))

# show growth rates correlate with expression response correlation
expression_response_growth <- inner_join(expression_response_correlation, growth_rates) %>%
  group_by(cell_id, drug) %>%
  top_n(-1, grmax) %>%
  ungroup() %>%
  mutate(cor_groups = cut(cor, c(-.4, -.3,-.2, -.1, 0, .1, .2)),
         cor_quant = cut(cor,
                         breaks = quantile(cor, probs = seq(0, 1, by = 0.2)),
                         labels = c("20","40","60","80","100"),
                         include.lowest=TRUE))

with(expression_response_growth,
     boxplot(grmax ~ cor_quant,
             xlab = 'Expression-response Correlation\n (Quantile)',
             ylab = 'Max Growth Inhibition\n (GRmax)',
             main = 'Drug Effects on Expression\n and Growth Rates'))
abline(h = 0, lty = 2, col = 'red')
dev.off()

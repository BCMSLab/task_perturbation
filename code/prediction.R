# load libraries ----
library(tidyverse)
library(clusterProfiler)

# ggplot theme ----
my_theme <- theme(panel.grid = element_blank(),
                  panel.border = element_rect(fill = NA, color = 'black', size = 1),
                  panel.background = element_blank(),
                  axis.title.x.top = element_blank(),
                  axis.text.x.top = element_blank(),
                  axis.title.y.right = element_blank(),
                  axis.text.y.right = element_blank(),
                  legend.background = element_blank(),
                  axis.ticks.length=unit(2, "mm"))

# load data ----
proteins <- read_rds('data_clean/breast_cells_proteins.rds')
genes <- read_rds('data_clean/breast_cells_genes.rds')
growth <- read_rds('data_clean/breast_cells_growth.rds')
hallmarks <- read_rds('data_clean/hallmarks_sets.rds')

common_cells <- intersect(colnames(genes), colnames(proteins))
common_symbols <- intersect(rownames(genes), rownames(proteins))

genes <- genes[common_symbols, common_cells]
proteins <- proteins[common_symbols, common_cells]

dim(genes); dim(proteins)

all(colnames(genes) == colnames(genes))
all(rownames(genes) == rownames(proteins))

# make the figure ----
gene_protein_corr <- melt(cor((genes), (proteins), use = 'complete')) %>%
  as_tibble() %>%
  mutate(group = ifelse(Var1 == Var2, 'Same cell', 'Other cells'))

gene_protein_plot <- gene_protein_corr %>%
  ggplot(aes(x = value, fill = group)) +
  geom_bar(stat="bin", aes(y=..density..), alpha = .5) +
  labs(x = 'Gene-protein correlation (PCC)',
       y = 'Percent of cell lines',
       fill = '') +
  scale_x_continuous(sec.axis = dup_axis(), limits = c(0, .4)) +
  scale_y_continuous(sec.axis = dup_axis()) +
  theme(legend.position = c(.35, .85)) +
  my_theme

gene_protein_hallmark <- map_df(hallmarks,
    function(x) {
      as_tibble(
        melt(
          cor(
            genes[intersect(rownames(genes), x),],
            proteins[intersect(rownames(proteins), x),],
            use = 'complete'
          )
        )
      )
    }, .id = 'hallmark') %>%
  filter(Var1 == Var2)

# enrichment of hallmakrs in breast cancer cells (genes)
term_gene <- melt(hallmarks) %>%
  as_tibble() %>%
  dplyr::select(term = L1, gene = value)

stats <- lapply(seq_len(ncol(genes)), function(i) genes[,i])
stats <- map(stats, sort, decreasing = TRUE)
names(stats) <- colnames(genes)

enrich_hallmark_genes <- map_df(stats, function(x) {
  enrich <- GSEA(x,
                 pAdjustMethod = 'fdr',
                 TERM2GENE = term_gene,
                 pvalueCutoff = 1,
                 maxGSSize = Inf)
  as_tibble(enrich)
}, .id = 'cell_id')

gene_protein_hallmark_plot <- gene_protein_hallmark %>%
  full_join(enrich_hallmark_genes,
            by = c('hallmark' = 'ID',
                   'Var1' = 'cell_id')) %>%
  dplyr::select(hallmark, Var1, value, NES) %>%
  mutate(hallmark = as.numeric(fct_reorder(hallmark, (value)))) %>%
  gather(type, newval, -hallmark, -Var1) %>%
  ggplot(aes(x = hallmark, y = newval, group = hallmark)) +
  geom_boxplot() +
  facet_wrap(~type, ncol = 1,
             scales = 'free_y',
             strip.position = "left", 
             labeller = as_labeller(c(value = 'Transcription-translation Correlation',
                                      NES = 'Enrichment Score'))) +
  labs(x = 'Cancer Hallmarks', y = '') +
  scale_x_continuous(sec.axis = dup_axis(),
                     breaks = 1:10) +
  scale_y_continuous(sec.axis = dup_axis()) +
  my_theme + 
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_text(size = 0))

plot_grid(
  plot_grid(NULL,
            gene_protein_plot,
            rel_heights = c(.7, 1),
            ncol = 1,
            scale = .9,
            labels = c('A', 'B'),
            label_fontface = 'plain'),
  gene_protein_hallmark_plot,
  nrow = 1,
  scale = c(1, .95),
  labels = c('', 'C'),
  label_fontface = 'plain'
) %>%
  ggsave(plot = .,
         filename = 'output/manuscript/figures/gene_protein.png',
         width = 5.6, height = 4.8)

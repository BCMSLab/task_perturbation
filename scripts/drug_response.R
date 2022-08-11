# Calculate the drug response in gene expression

# load libraries
library(tidyverse)
library(slinky)
library(org.Hs.eg.db)
library(SummarizedExperiment)

# download.file('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_inst_info.txt.gz',
#               destfile = 'data/GSE92742_Broad_LINCS_inst_info.txt.gz')
# aria2c -x 8 -s 8 https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl//GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx.gz

# link to gctx files and key
user_key <- httr::content(httr::GET("https://api.clue.io/temp_api_key"), as = "parsed")$user_key
gctx_file <- 'data/GSE92742_Broad_LINCS_Level3_INF_mlr12k_n1319138x12328.gctx'
gctx_info <- 'data/GSE92742_Broad_LINCS_inst_info.txt.gz'

# load gctx file and metadata
sl <- Slinky(user_key, gctx_file, gctx_info)
se <- read_rds('output/gene_expression.rds')

# for(cell in unique(se$cell_id)) {
#   lapply(unique(se$pert_iname),
#            function(drug) {
#              n <- sum(se$cell_id == cell & se$pert_iname == drug)
#              if (n > 5) {
#                ks_vecs <- diffexp(
#                  sl,
#                  treat = se[, se$cell_id == cell & se$pert_iname == drug],
#                  control = se[, se$cell_id == cell & se$pert_iname == 'DMSO'],
#                  method = 'ks',
#                  verbose = TRUE
#                )
#                
#                df <- melt(ks_vecs)
#                names(df) <- c('gene', drug)
#                df$cell_id <- cell
#                
#                write_csv(df, file.path('scores', paste0(cell, '_', drug, '.csv')))
#              }
#            })
# }

fls <- list.files('data/scores', full.names = TRUE)
names(fls) <- str_split(fls, '\\.|/', simplify = TRUE)[, 3]

drug_response <- imap(fls, function(x, .y) {
  df <- read_csv(x)[, 2]
  names(df) <- .y
  df
})

drug_response <- do.call(cbind, drug_response) %>% as_tibble()

entrez <- read_csv(fls[1])[,1,drop=TRUE]

id_symbol <- AnnotationDbi::select(
  org.Hs.eg.db,
  as.character(entrez),
  'SYMBOL',
  'ENTREZID')
m <- as.matrix(drug_response)
rownames(m) <- id_symbol$SYMBOL

pheno_data <- tibble(file = colnames(m)) %>%
  separate(file, into = c('cell_id', 'drug')) %>%
  as.data.frame()

se_response <- SummarizedExperiment(assays = list(response = m),
                                    colData = pheno_data)

write_rds(se_response, 'output/drug_response.rds')

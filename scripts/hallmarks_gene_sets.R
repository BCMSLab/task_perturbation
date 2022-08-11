# Download the cancer hallmarks gene sets

# load libraries
library(tidyverse)
library(reshape2)
library(readxl)
library(org.Hs.eg.db)

# hallmarks
# source: CHG: A Systematically Integrated Database of Cancer Hallmark Genes
# source: http://bio-bigdata.hrbmu.edu.cn/CHG/index.html
url <- 'http://bio-bigdata.hrbmu.edu.cn/CHG/download/Supplementary%20Table%202%20Genesets%20of%2010%20hallmarker.xlsx'

download.file(url, 'data/hallmarks_gene_sets.xlsx')
symbols <- keys(org.Hs.eg.db, 'SYMBOL')

# read file
hallmarks <- read_excel('data/hallmarks_gene_sets.xlsx', skip = 2)

# transform symbols and subset to valid
hallmarks <- as.list(hallmarks) %>%
  map(function(x) {
    l <- x[!is.na(x)]
    l <- toupper(l)
    intersect(l, symbols)
  })

# tidy
hallmarks <- melt(hallmarks) %>% as_tibble()
names(hallmarks) <- c('symbol', 'hallmark')

write_csv(hallmarks, 'output/hallmarks_gene_sets.csv')

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(readr)
library(stringr)
library(FactoMineR)
library(factoextra)
library(scrime)
library(Hmisc)
library(purrr)
library(uwot)

# RF pre-processed dataset
snp_tbl = readr::read_csv(file = 'data/snp_tbl.csv')
metadata = snp_tbl %>% select(c(1,2,3,4)) # keep just the metadata

snp_rf_data = readr::read_csv(file = 'data/rf_preprocessed_snp_data.csv', col_names = TRUE)
all(snp_rf_data$poparea == metadata$poparea)

X = snp_rf_data %>% select(-c(1,2)) # remove poparea and sample id

# data density
# X %>% unname() %>% unlist() %>% as.vector() %>% density(adjust = 0.3) %>% plot()

# PCA
set.seed(42)
pca_res = FactoMineR::PCA(X, scale.unit = FALSE, graph = FALSE) # fast

## 4 colors for geography and 12 colors for population area with 3:1 distinction
## ratio (every 3 population colors are alike, corresponding to the geography)
geo_col = RColorBrewer::brewer.pal(n = 9, name = 'Set1')[c(1,2,3,9)]
pop_col = c('#8c1011', '#e41a1c', '#ef7172',
            '#204a6c', '#377eb8', '#89b7dc',
            '#295c27', '#4DAF4A', '#9cd49a',
            '#414141', '#999999', '#cacaca')

factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = metadata$geography %>% factor(levels = c('RA', 'VE', 'MD', 'NN')),
  addEllipses = TRUE,
  legend.title = "Geography", pointshape = 20, palette = geo_col,
  title = "PCA")
ggsave('fig/rf_pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = metadata$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_pca_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_pca_pop.pdf', width = 7, height = 5)

# UMAP
neighbors = c(2,4,6,8,10,13,15,20,24)
metrics = c('euclidean', 'manhattan', 'cosine')

for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_res = uwot::umap(X = X, n_neighbors = nn, metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$geography %>%
          factor(levels = c('RA', 'VE', 'MD', 'NN'))) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col) +
      labs(title = paste0("UMAP (RF-preprocessed SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_umap_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (RF-preprocessed SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_umap_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (RF-preprocessed SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/rf_umap_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}

# read RF-reduced dataset
snp_rf_data_reduced = readr::read_csv(file = 'data/rf_final_330snp.csv', col_names = TRUE)
all(snp_rf_data_reduced$poparea == metadata$poparea)

X_red = snp_rf_data_reduced %>% select(-c(1,2)) # remove poparea and sample id

set.seed(42)
pca_res2 = FactoMineR::PCA(X_red, scale.unit = FALSE, graph = FALSE) # fast

# PCA
factoextra::fviz_pca_ind(pca_res2, geom.ind = "point",
  col.ind = metadata$geography %>% factor(levels = c('RA', 'VE', 'MD', 'NN')),
  addEllipses = TRUE, palette = geo_col,
  legend.title = "Geography", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_red_pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res2, geom.ind = "point",
  col.ind = metadata$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_red_pca_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res2, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_red_pca_pop.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res2, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = FALSE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_red_pca_pop_no_ellipses.pdf', width = 7, height = 5)

# UMAP
for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_res = uwot::umap(X = X_red, n_neighbors = nn, metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$geography %>%
          factor(levels = c('RA', 'VE', 'MD', 'NN'))) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col) +
      labs(title = paste0("UMAP (RF-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_red_umap_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (RF-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_red_umap_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (RF-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/rf_red_umap_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}

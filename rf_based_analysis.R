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

########################################################################
# RF pre-processed dataset (~10000 features, for Geography prediction) #
########################################################################
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

##################################################################
# RF rf-reduced dataset (330 features, for Geography prediction) #
##################################################################
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

############################################################
# RF rf-reduced dataset (371 features, for Sex prediction) #
############################################################
snp_rf_rex_reduced = readr::read_csv(file = 'data/rf_sex_371snp_final.csv', col_names = TRUE)
all(snp_rf_rex_reduced$poparea == metadata$poparea)
all(snp_rf_rex_reduced$sex == metadata$sex)
all(snp_rf_rex_reduced$geography == metadata$geography)

X_red3 = snp_rf_rex_reduced %>% select(-c(1,2,3,4)) # remove sample id, poparea, sex and geography

set.seed(42)
pca_res3 = FactoMineR::PCA(X_red3, scale.unit = FALSE, graph = FALSE) # fast

# PCA
factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$geography %>% factor(levels = c('RA', 'VE', 'MD', 'NN')),
  addEllipses = TRUE, palette = geo_col,
  legend.title = "Geography", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_pop.pdf', width = 7, height = 5)

# UMAP
neighbors = c(2,8,13,15,20)
metrics = c('euclidean', 'manhattan', 'cosine')

for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_res = uwot::umap(X = X_red3, n_neighbors = nn, metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$geography %>%
          factor(levels = c('RA', 'VE', 'MD', 'NN'))) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col) +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}
############################################################
# RF rf-reduced dataset (371 features, for Sex prediction) #
############################################################
snp_rf_rex_reduced = readr::read_csv(file = 'data/rf_sex_371snp_final.csv', col_names = TRUE)
all(snp_rf_rex_reduced$poparea == metadata$poparea)
all(snp_rf_rex_reduced$sex == metadata$sex)
all(snp_rf_rex_reduced$geography == metadata$geography)

X_red3 = snp_rf_rex_reduced %>% select(-c(1,2,3,4)) # remove sample id, poparea, sex and geography

set.seed(42)
pca_res3 = FactoMineR::PCA(X_red3, scale.unit = FALSE, graph = FALSE) # fast

# PCA
factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$geography %>% factor(levels = c('RA', 'VE', 'MD', 'NN')),
  addEllipses = TRUE, palette = geo_col,
  legend.title = "Geography", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res3, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_sex_red_pca_pop.pdf', width = 7, height = 5)

# UMAP
neighbors = c(2,8,13,15,20)
metrics = c('euclidean', 'manhattan', 'cosine')

for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_res = uwot::umap(X = X_red3, n_neighbors = nn, metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$geography %>%
          factor(levels = c('RA', 'VE', 'MD', 'NN'))) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col) +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (RF-sex-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/rf_sex_red_umap_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}

################################################################
# RF rf-reduced dataset (267 features, for Poparea prediction) #
################################################################
snp_rf_poparea_reduced = readr::read_csv(file = 'data/rf_poparea_267snp_final.csv', col_names = TRUE)
all(snp_rf_poparea_reduced$poparea == metadata$poparea)
all(snp_rf_poparea_reduced$sex == metadata$sex)
all(snp_rf_poparea_reduced$geography == metadata$geography)

X_red4 = snp_rf_poparea_reduced %>% select(-c(1,2,3,4)) # remove sample id, poparea, sex and geography

set.seed(42)
pca_res4 = FactoMineR::PCA(X_red4, scale.unit = FALSE, graph = FALSE) # fast

# PCA
factoextra::fviz_pca_ind(pca_res4, geom.ind = "point",
  col.ind = metadata$geography %>% factor(levels = c('RA', 'VE', 'MD', 'NN')),
  addEllipses = TRUE, palette = geo_col,
  legend.title = "Geography", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_poparea_red_pca_geo.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res4, geom.ind = "point",
  col.ind = metadata$sex %>% as.factor(), addEllipses = TRUE,
  legend.title = "Sex", pointshape = 20,
  title = "PCA")
ggsave('fig/rf_poparea_red_pca_sex.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res4, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = TRUE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_poparea_red_pca_pop.pdf', width = 7, height = 5)

factoextra::fviz_pca_ind(pca_res4, geom.ind = "point",
  col.ind = metadata$poparea %>% as.factor(), addEllipses = FALSE,
  legend.title = "Pop", pointshape = 20, palette = pop_col,
  title = "PCA")
ggsave('fig/rf_poparea_red_pca_pop_no_ellipses.pdf', width = 7, height = 5)

# UMAP
neighbors = c(2,8,13,15,20)
metrics = c('euclidean', 'manhattan', 'cosine')

for (nn in neighbors) {
  for (met in metrics) {
    print(paste0('Neighbors: ', nn, ", Metric: ", met))

    set.seed(42)
    umap_res = uwot::umap(X = X_red4, n_neighbors = nn, metric = met, verbose = TRUE)

    # UMAP figure (Geography coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$geography %>%
          factor(levels = c('RA', 'VE', 'MD', 'NN'))) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = geo_col) +
      labs(title = paste0("UMAP (RF-poparea-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Geography",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_poparea_red_umap_', nn, '_', met,'_geo.pdf'), width = 7, height = 5)

    # UMAP figure (Sex coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$sex %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      labs(title = paste0("UMAP (RF-poparea-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Sex",
        label.theme = element_text(size = 12),
        override.aes = list(shape = 19, size = 12)))
    ggsave(paste0('fig/rf_poparea_red_umap_', nn, '_', met,'_sex.pdf'), width = 7, height = 5)

    # UMAP figure (Population coloring)
    umap_res %>%
      `colnames<-` (c("X", "Y")) %>%
      tibble::as_tibble() %>%
      tibble::add_column(cluster = metadata$poparea %>% as.factor()) %>%
      ggplot(aes(x = X, y = Y, color = cluster)) +
      geom_point() +
      scale_color_manual(values = pop_col) +
      labs(title = paste0("UMAP (RF-poparea-reduced SNP data, ", nn, " Neighbors, ", met, ")")) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      guides(color = guide_legend(title = "Pop",
        label.theme = element_text(size = 10),
        override.aes = list(size = 3)))
    ggsave(paste0('fig/rf_poparea_red_umap_', nn, '_', met,'_pop.pdf'), width = 7, height = 5)
  }
}

# Venn Digramm comparing SNP features of reduced RF-datasets
set1_col = RColorBrewer::brewer.pal(n = 3, name = 'Set1')
geo_features = colnames(X_red)
sex_features = colnames(X_red3)
poparea_features = colnames(X_red4)

# supress logger output
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

VennDiagram::venn.diagram(x = list(geo_features, sex_features, poparea_features),
  category.names = c("Geography" , "Sex" , "Population Area"),
  main = 'Common SNP features',
  filename = 'fig/venn.png', output = TRUE, imagetype = "png",
  lty = 'blank', fill = set1_col, cex = 2, margin = 0.1, cat.cex = 1.6)

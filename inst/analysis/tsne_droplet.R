library(Matrix)
library(Rtsne)
library(ggplot2)
library(cowplot)
set.seed(1)
load("results.RData")
load("~/git/fastTopics-experiments/data/droplet.RData")
V1 <- droplets_res_list$scGBM[["10_factors"]]$V
V2 <- droplets_res_list$glmpca[["10_factors"]]$V
V3 <- droplets_res_list$fastglmpca_28_cores$PCs_10$V

tissue_colors <- c("royalblue",   # basal
                   "firebrick",   # ciliated
                   "forestgreen", # club
                   "gold",        # goblet
                   "darkmagenta", # ionocyte
                   "darkorange",  # neuroendocrine
                   "skyblue")     # tuft

# Compute 2-d t-SNEs from GLM-PCA fits.
run_tsne <- function (V) {
  out <- Rtsne(V,dims = 2,pca = FALSE,normalize = FALSE,eeta = 200,
               perplexity = 100,theta = 0.1,max_iter = 1000,
               check_duplicates = FALSE,verbose = TRUE,
               num_threads = 4)
  colnames(out$Y) <- c("d1","d2")
  return(out$Y)
}
Y1 <- run_tsne(V1)
Y2 <- run_tsne(V2)
Y3 <- run_tsne(V3)

# Create the t-SNE plots.
pdat1 <- cbind(samples,Y1)
pdat2 <- cbind(samples,Y2)
pdat3 <- cbind(samples,Y3)
p1 <- ggplot(pdat1,aes(x = -d1,y = d2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("scGBM") +
  theme_cowplot(font_size = 9)
p2 <- ggplot(pdat2,aes(x = d2,y = -d1,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("glmpca") +
  theme_cowplot(font_size = 9)
p3 <- ggplot(pdat3,aes(x = d1,y = d2,color = tissue)) +
  geom_point(size = 1) +
  scale_color_manual(values = tissue_colors) +
  ggtitle("fastglmpca") +
  theme_cowplot(font_size = 9)
plot_grid(p1,p2,p3,nrow = 1,ncol = 3)
ggsave("tsne.eps",
       plot_grid(p1,p2,p3,nrow = 1,ncol = 3),
       height = 2,width = 9.25)

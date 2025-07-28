# Title: R script for the publication (DOI:)
# Author: Gema Rodríguez-Caballero
# Email: a62rocag@uco.es
# License: MIT
# Date: 2025-07-28
  

# --- Load data ---
ASV <- read.table(file = "ASVs_TotalAbund.csv", header = TRUE, sep = ";", dec = ",")
row.names(ASV) <- ASV$Sample

# --- Rarefaction ---
library(vegan)

set.seed(4)
ASV_rar <- rrarefy(ASV[, -c(1:2)], min(rowSums(ASV[, -c(1:2)])))

# --- Rarefaction Curves ---
species_counts <- specnumber(ASV_rar)
raremax <- min(rowSums(ASV_rar))
S_rare <- rarefy(ASV_rar, raremax)

# Basic rarefaction plot
plot(species_counts, S_rare,
     xlab = "Observed No. of Species",
     ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(ASV_rar, step = 20, sample = raremax, col = "blue", cex = 0.6)


# --- PERMANOVA ---
set.seed(15)
perMANOVA <- adonis2(ASV_rar ~ PlasticType, data = ASV, perm = 999, distance = "bray")
print(perMANOVA)

# Pairwise PERMANOVA
library(pairwiseAdonis)
set.seed(94)
PW_adonis <- pairwise.adonis(ASV_rar, ASV$PlasticType, p.adjust.m = "BH", sim.method = "bray")
print(PW_adonis)

# --- NMDS ---
set.seed(14)
NMDS_result <- metaMDS(ASV_rar, k = 2, wascores = TRUE, distance = "bray")
stressplot(NMDS_result)
plot(NMDS_result)

# Prepare NMDS plot data
NMDS_data <- as.data.frame(NMDS_result$points)
NMDS_data$PlasticType <- as.factor(ASV$PlasticType)
NMDS_data$Sample <- ASV$Sample

# --- NMDS Plot with Ellipses ---
library(ggplot2)
library(ggrepel)

# NMDS base plot
plot_NMDS <- ggplot(data = NMDS_data, aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = PlasticType), size = 1.5) +  
  geom_text_repel(aes(label = Sample)) +
  scale_size(guide = "none") +
  theme_classic() +
  theme(panel.border = element_rect(linetype = "solid", colour = "black", fill = NA)) +
  xlab("NMDS1") +
  ylab("NMDS2")

# Custom ellipse function (from https://stackoverflow.com/a/13801089 by Didzis Elferts https://www.delferts.lv/en/)
veganCovEllipse <- function(cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- seq(0, 2 * pi, length.out = npoints)
  circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(circle %*% chol(cov)))
}

# Calculate ellipses (requires `ordiellipse` result)
ellip <- ordiellipse(NMDS_result, ASV$PlasticType, display = "sites", kind = "sd", label = FALSE)
df_ell <- data.frame()
for (g in levels(NMDS_data$PlasticType)) {
  ell <- with(NMDS_data[NMDS_data$PlasticType == g, ],
              veganCovEllipse(ellip[[g]]$cov, ellip[[g]]$center, ellip[[g]]$scale))
  df_ell <- rbind(df_ell, cbind(as.data.frame(ell), PlasticType = g))
}
colnames(df_ell)[1:2] <- c("NMDS1", "NMDS2")

# Add ellipses to NMDS plot
plot_NMDS <- plot_NMDS +
  geom_path(data = df_ell, aes(x = NMDS1, y = NMDS2, colour = PlasticType),
            size = 1, linetype = 2)

# --- Hierarchical Clustering ---
ASV2 <- ASV[, -c(1:2)]

# Bray-Curtis distance and Ward.D2 clustering
dist_BC <- vegdist(ASV_rar, method = "bray")
hclust_result <- hclust(dist_BC, method = "ward.D2")
dendrogram <- as.dendrogram(hclust_result)

# Plot dendrogram
maxHeight <- max(hclust_result$height)

plot(dendrogram, cex = 0.5, ylim = c(-maxHeight * 0.2, maxHeight), horiz = FALSE, type = "r")
rect.hclust(hclust_result, k = 3)


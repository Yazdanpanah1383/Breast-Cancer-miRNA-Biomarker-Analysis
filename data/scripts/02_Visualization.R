#############################################################
# Breast Cancer miRNA Biomarker Analysis
# Step 2: Data Visualization
#############################################################

library(ggplot2)
library(pheatmap)

#############################################################
# Volcano Plot
#############################################################

# Prepare adjusted P-values for plotting

res$plotP <- pmin(-log10(res$adj.P.Val + 1e-300), 50)
sig$plotP <- pmin(-log10(sig$adj.P.Val + 1e-300), 50)

volcano_plot <- ggplot() +
  
  geom_point(
    data = res,
    aes(x = logFC, y = plotP),
    color = "grey70",
    alpha = 0.6,
    size = 2
  ) +
  
  geom_point(
    data = sig,
    aes(x = logFC, y = plotP),
    color = "red",
    size = 2.5
  ) +
  
  geom_vline(
    xintercept = c(-1,1),
    linetype = "dashed",
    color = "blue"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "blue"
  ) +
  
  labs(
    title = "Volcano Plot",
    x = expression(log[2]("Fold Change")),
    y = expression(-log[10]("Adjusted P-value"))
  ) +
  
  coord_cartesian(ylim = c(0,50)) +
  
  theme_classic(base_size = 15)

volcano_plot

#############################################################
# Heatmap
#############################################################

heat_data <- expr[top24$ID_REF, ]

pheatmap(
  heat_data,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE
)

#############################################################
# PCA Plot
#############################################################

expr[!is.finite(expr)] <- 0

pca <- prcomp(
  t(expr),
  scale. = TRUE
)

group <- rep(NA, ncol(expr))

group[bc1_start:bc1_end] <- "Cancer"
group[nc_start:nc_end] <- "Normal"
group[ben_start:ben_end] <- "Benign"

pca_data <- data.frame(
  
  PC1 = pca$x[,1],
  
  PC2 = pca$x[,2],
  
  Group = group
  
)

pca_plot <- ggplot(
  pca_data,
  aes(PC1, PC2, color = Group)
) +
  
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  stat_ellipse(level = 0.95) +
  
  theme_classic(base_size = 15) +
  
  labs(
    title = "PCA Plot",
    x = "PC1",
    y = "PC2"
  )

pca_plot

#############################################################
# Boxplot
#############################################################

box_df <- data.frame(
  
  Expression = as.numeric(expr[top5$ID_REF[1], ]),
  
  Group = factor(group)
  
)

boxplot1 <- ggplot(
  
  box_df,
  
  aes(Group, Expression, fill = Group)
  
) +
  
  geom_boxplot(alpha = 0.7) +
  
  geom_jitter(width = 0.2) +
  
  theme_classic(base_size = 15) +
  
  labs(
    
    title = paste("Expression of", top5$ID_REF[1]),
    
    x = "",
    
    y = "Expression"
    
  )

boxplot1

#############################################################
# Save Figures
#############################################################

ggsave("Volcano.png", volcano_plot,
       width = 8,
       height = 6,
       dpi = 300)

ggsave("PCA.png", pca_plot,
       width = 7,
       height = 6,
       dpi = 300)

ggsave("Boxplot.png", boxplot1,
       width = 6,
       height = 5,
       dpi = 300)

#############################################################
# End of Script
#############################################################

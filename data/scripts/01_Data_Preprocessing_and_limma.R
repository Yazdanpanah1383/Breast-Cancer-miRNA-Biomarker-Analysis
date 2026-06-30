#############################################################
# Breast Cancer miRNA Biomarker Analysis
# Step 1: Data Preprocessing & Differential Expression
#############################################################

# Clear workspace 
rm(list = ls())

# ===========================
# Load Packages
# ===========================

library(limma)

# ===========================
# Read Data
# ===========================

# Read data (CSV file)
# data <- read.csv("YourData.csv", check.names = FALSE)

# Skip this step if the dataset has already been loaded into the R environment.
# ===========================
# Data Structure
# ===========================

class(data)

str(data)

# # Convert columns to numeric

data[,-1] <- lapply(data[,-1], as.numeric)

expr <- as.matrix(data[,-1])

rownames(expr) <- data$ID_REF

mode(expr) <- "numeric"

# ===========================
# Group Information
# ===========================

group <- factor(c(
  
  rep("Cancer", length(bc1_start:bc1_end)),
  
  rep("Normal", length(nc_start:nc_end)),
  
  rep("Benign", length(ben_start:ben_end))
  
))

design <- model.matrix(~0 + group)

colnames(design) <- levels(group)

design

# ===========================
# Contrast Matrix
# ===========================

contrast.matrix <- makeContrasts(
  
  Cancer_vs_Normal = Cancer - Normal,
  
  levels = design
  
)

contrast.matrix

# ===========================
# limma Analysis
# ===========================

fit <- lmFit(expr, design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

# ===========================
# Differential Expression
# ===========================

res <- topTable(
  
  fit2,
  
  coef = "Cancer_vs_Normal",
  
  number = Inf
  
)

# Add miRNA identifiers

res$ID_REF <- data$ID_REF

# Display the results

head(res)

head(res[,c("ID_REF","logFC","adj.P.Val")])

# ===========================
# Significant miRNAs
# ===========================

sig <- subset(
  
  res,
  
  adj.P.Val < 0.001 &
    
    abs(logFC) > 3
  
)

cat("Number of significant miRNAs:", nrow(sig), "\n")

# ===========================
# Top24
# ===========================

top24 <- sig

write.csv(
  
  top24,
  
  "Top24_miRNAs.csv",
  
  row.names = FALSE
  
)

# ===========================
# Top5
# ===========================

top5 <- head(
  
  sig[order(sig$adj.P.Val), ],
  
  5
  
)

top5

#############################################################
# End of File 01
#############################################################

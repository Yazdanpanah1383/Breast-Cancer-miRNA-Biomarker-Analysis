#############################################################
# Breast Cancer miRNA Biomarker Analysis
# Step 3: ROC Curve Analysis
#############################################################

library(pROC)

#############################################################
# Prepare ROC Data
#############################################################

# Cancer = 1
# Normal = 0
# Benign samples are excluded

roc_group <- ifelse(group == "Cancer", 1,
                    ifelse(group == "Normal", 0, NA))

keep <- !is.na(roc_group)

#############################################################
# Top 5 miRNAs
#############################################################

top5 <- head(sig[order(sig$adj.P.Val), ], 5)

#############################################################
# Individual ROC Curves
#############################################################

auc_results <- data.frame(
  miRNA = character(),
  AUC = numeric()
)

for(i in 1:nrow(top5)){
  
  mir <- top5$ID_REF[i]
  
  expression <- as.numeric(expr[mir, keep])
  
  roc_obj <- roc(
    response = roc_group[keep],
    predictor = expression,
    quiet = TRUE
  )
  
  auc_results <- rbind(
    auc_results,
    data.frame(
      miRNA = mir,
      AUC = as.numeric(auc(roc_obj))
    )
  )
  
  png(
    paste0("ROC_", mir, ".png"),
    width = 1800,
    height = 1800,
    res = 300
  )
  
  plot(
    roc_obj,
    col = "#D55E00",
    lwd = 3,
    main = paste("ROC Curve -", mir)
  )
  
  legend(
    "bottomright",
    legend = paste("AUC =", round(auc(roc_obj),3)),
    bty = "n"
  )
  
  dev.off()
  
}

#############################################################
# Combined ROC Curve
#############################################################

combined_score <- colMeans(expr[top5$ID_REF, keep])

roc_combined <- roc(
  response = roc_group[keep],
  predictor = combined_score,
  quiet = TRUE
)

png(
  "ROC_Top5_Combined.png",
  width = 1800,
  height = 1800,
  res = 300
)

plot(
  roc_combined,
  col = "blue",
  lwd = 3,
  main = "Combined ROC Curve"
)

legend(
  "bottomright",
  legend = paste(
    "AUC =",
    round(auc(roc_combined),3)
  ),
  bty = "n"
)

dev.off()

#############################################################
# Save AUC Table
#############################################################

write.csv(
  auc_results,
  "Top5_AUC_Results.csv",
  row.names = FALSE
)

#############################################################
# Print Combined AUC
#############################################################

cat(
  "Combined AUC:",
  round(auc(roc_combined),3),
  "\n"
)

#############################################################
# End of Script
#############################################################
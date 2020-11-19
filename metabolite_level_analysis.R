# This script depicts heatmap plotting of the metabolite level fold change 

load("~/peak_filtered_cleanup.Rda")

# Fold change calculation
peak_filter3$col0_fold_change <- peak_filter3$col_treated_mean/peak_filter3$col_ctrl_mean
peak_filter3$idm1_fold_change <- peak_filter3$idm_treated_mean/peak_filter3$idm_ctrl_mean
peak_filter3$clf28_fold_change <- peak_filter3$clf_treated_mean/peak_filter3$clf_ctrl_mean
peak_filter3_sort <- peak_filter3[order(-peak_filter3$clf28_fold_change),]

fold_change_plot <- subset(peak_filter3_sort, select = c(m_z, col0_fold_change, idm1_fold_change, clf28_fold_change))

# Heatmap plotting
tiff(filename="clf_fold_change_heatmap.tiff", height = 100, width = 5, units = 'in', res=600)
col <- colorRampPalette(c("orange","white","darkgreen"))(200)
heatmap.2(as.matrix(fold_change_plot[,2:4]), Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "none", col = col, density.info = "none", trace = "column", tracecol = "black", key = 0)
dev.off()

# Selected inspection of benzoic acid and camalexin
selected_metabolites <- fold_change_plot[c("123.09947148","201.14684081"),]
heatmap.2(as.matrix(selected_metabolites[,2:4]), Rowv = FALSE, Colv = FALSE, dendrogram = "none", scale = "none", col = col, density.info = "none", trace = "column", tracecol = "black", key = 0)

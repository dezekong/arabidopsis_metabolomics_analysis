# This script depicts the data clean up process for 'data_table' obtained from 'peak_extraction_from_mzxml_files.R'

load("~/data_table.Rda")

# Re-order the columns
col_order <- c("peak_names","mz_med","mz_min","mz_max","rt_med","rt_min","rt_max","n_peaks_ctrl_treated",
               "100nM_cam_std","col_ctrl1","col_ctrl2","col_ctrl3","col_ctrl4","col_treated1","col_treated2","col_treated3",
               "clf_ctrl1","clf_ctrl2","clf_ctrl3","clf_ctrl4","clf_treated1","clf_treated2","clf_treated3","clf_treated4",
               "idm_ctrl1","idm_ctrl2","idm_ctrl3","idm_ctrl4","idm_treated1","idm_treated2","idm_treated3","idm_treated4",
               "100nM_cam_std_int","col_ctrl1_int","col_ctrl2_int","col_ctrl3_int","col_ctrl4_int","col_treated1_int","col_treated2_int","col_treated3_int",
               "clf_ctrl1_int","clf_ctrl2_int","clf_ctrl3_int","clf_ctrl4_int","clf_treated1_int","clf_treated2_int","clf_treated3_int","clf_treated4_int",
               "idm_ctrl1_int","idm_ctrl2_int","idm_ctrl3_int","idm_ctrl4_int","idm_treated1_int","idm_treated2_int","idm_treated3_int","idm_treated4_int")
data_table <- data_table[, col_order]

# Calculate the total number of peaks detected for each metabolite across all samples - the n_peak parameter from original dataset is not accurate
data_table$n_peaks_ctrl_treated <- rowSums(data_table[,10:32])

# Check the data collected for camalexin, obtain retention time info from authentic std
# Confirm camalexin detection in all samples, remove samples with failed detection (if any)
camalexin_candidates <- subset(data_table, mz_med < 202 & mz_med > 201)

# Sort the data table based on n_peaks_ctrl_treated in descending order
data_table_sort <- data_table[order(-data_table$n_peaks_ctrl_treated),]

# Select data rows with more than 21 peaks, and no data at "idm_treated1", due to loss of sample during measurement
data_table_22_peaks <- data_table_sort[data_table_sort[,"n_peaks_ctrl_treated"] > 21,]
data_table_22_peaks <- data_table_22_peaks[data_table_22_peaks[,"idm_treated1"] == 0,]

# Calculate peak intensity mean among triplicates - data quality monitored by peak FT028939 (camalexin)
data_table_22_peaks$col_ctrl_mean <- rowMeans(data_table_22_peaks[,c('col_ctrl1_int','col_ctrl3_int','col_ctrl4_int')], na.rm=TRUE)
data_table_22_peaks$col_treated_mean <- rowMeans(data_table_22_peaks[,c('col_treated1_int', 'col_treated2_int','col_treated3_int')], na.rm=TRUE)
data_table_22_peaks$clf_ctrl_mean <- rowMeans(data_table_22_peaks[,c('clf_ctrl1_int', 'clf_ctrl2_int','clf_ctrl4_int')], na.rm=TRUE)
data_table_22_peaks$clf_treated_mean <- rowMeans(data_table_22_peaks[,c('clf_treated2_int', 'clf_treated3_int','clf_treated4_int')], na.rm=TRUE)
data_table_22_peaks$idm_ctrl_mean <- rowMeans(data_table_22_peaks[,c('idm_ctrl1_int', 'idm_ctrl2_int','idm_ctrl4_int')], na.rm=TRUE)
data_table_22_peaks$idm_treated_mean <- rowMeans(data_table_22_peaks[,c('idm_treated2_int', 'idm_treated3_int','idm_treated4_int')], na.rm=TRUE)

# Apply filter to dataset: clf_treated_mean > clf_ctrl_mean
data_table_peak_filter3 <- data_table_22_peaks[data_table_22_peaks[,"clf_treated_mean"] > data_table_22_peaks[,"clf_ctrl_mean"],]

# Column-wise clean-up, rename row with m/z
row.names(data_table_peak_filter3) <- round(data_table_peak_filter3[,"mz_med"], digits = 8)
peak_filter3 <- subset(data_table_peak_filter3, select = c(col_ctrl1_int, col_ctrl3_int, col_ctrl4_int, col_ctrl_mean,
                                                           col_treated1_int, col_treated2_int, col_treated3_int, col_treated_mean,
                                                           clf_ctrl1_int, clf_ctrl2_int, clf_ctrl4_int, clf_ctrl_mean,
                                                           clf_treated2_int, clf_treated3_int, clf_treated4_int, clf_treated_mean,
                                                           idm_ctrl1_int, idm_ctrl2_int, idm_ctrl4_int, idm_ctrl_mean,
                                                           idm_treated2_int, idm_treated3_int, idm_treated4_int, idm_treated_mean))
peak_filter3$m_z <- round(data_table_peak_filter3[,"mz_med"], digits = 8)

save(peak_filter3,file="peak_filtered_cleanup.Rda")
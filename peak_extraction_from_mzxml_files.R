# This script depicts the analysis of Arabidopsis thaliana metabolomics data collected from Agilent qTOF 6545
# The Aligent .d files were centroided and converted to .mzXML files via MSConvert prior to analysis

library(xcms)
library(faahKO)
library(multtest)

# Collect detected .mzXML files under the defined directory
mzmlfiles <- list.files(full.names = TRUE, recursive = TRUE)

# Read the centroided data
msdata <- readMSData(mzmlfiles, msLevel = 1, mode = "onDisk")

# Pick the peaks - parameters optimized by the instrument specs & targeted peak size
cwp <- CentWaveParam(ppm = 20, peakwidth = c(5,30), integrate = 2, mzdiff = -0.001)
peaklist <- findChromPeaks(msdata, param = cwp)

# Group the peaks
pdp <- PeakDensityParam(sampleGroups = fileNames(peaklist), bw = 5, minFraction = 1, binSize = 0.001)
peaklist <- groupChromPeaks(peaklist, pdp, msLevel = 1L, add = FALSE)

# Alignment, retention time correction - for peak visualization
peaklist_rc <- peaklist
pgp <- PeakGroupsParam(minFraction = 1, smooth = "loess", span = 0.6, family = "gaussian")
peaklist_rc <- adjustRtime(peaklist_rc, param = pgp, msLevel = 1)
plotAdjustedRtime(peaklist_rc)

# Obtain correspondence & intensity results
peak_correspondence <- featureDefinitions(peaklist)
peak_intensity <- featureValues(peaklist, value = "into")

# Generate data table
data_table <- merge(peak_correspondence, peak_intensity, by=0, all=TRUE)

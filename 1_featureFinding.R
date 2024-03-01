library(patRoon)
library(xcms)
library(tidyverse)

options(patRoon.MP.maxProcs = 4) # Limit max processes in patRoon
options(future.globals.maxSize = 1000000000000) # increases max memory for MS Spectra generation (reporting)
register(SerialParam()) # Limit processes to prevent XCMS crashing (RT alignment & grouping)

# -------------------------
# Convert raw files to mzML
# -------------------------

convertMSFiles("RAW", "mzML", dirs = TRUE, from = "thermo",
               centroid = "vendor", filters = "precursorRecalculation",
               overWrite = TRUE)

# -------------------------
# Load analysis information
# -------------------------

# Manually created list from multiple paths
anaInfo <- read.csv("anaInfo.csv", check.names = FALSE)


# -------------------------
# Find features from all samples
# -------------------------

# Peak picking
time.xcmsPick <- system.time({
  fList <- findFeatures(
    anaInfo,
    "xcms3",
    param = xcms::CentWaveParam(
      ppm = 5,
      peakwidth = c(10, 60), # half avg peak width - 2x avg peak width
      snthresh = 3, # 10 last successful test
      prefilter = c(3, 100),
      mzCenterFun = "wMean", # from wMeanApex3
      integrate = 1L,
      mzdiff = 0.005, # minimum difference in m/z dimension required for peaks with overlapping retention times
      fitgauss = TRUE, # normally false
      noise = 1000
    ),
    verbose = FALSE
  )
})

saveRDS(fList, "fList-original.rds")
fList <- readRDS("fList-original.rds")

# Group features and align retention time (114965)
time.xcmsGroup <- system.time({
  fGroups <- groupFeatures(
    fList,
    "xcms3",
    rtalign = TRUE,
    loadRawData = TRUE,
    groupParam = xcms::PeakDensityParam(
      sampleGroups = anaInfo$group,
      minFraction = 0,
      minSamples = 1,
      bw = 15, # cranked from 10 due to late eluting big peaks
      binSize = 0.01 # corrected for misaligned m/z in features
    ),
    retAlignParam = xcms::ObiwarpParam(
      center = 5,
      response = 1,
      gapInit = 0.3, #0.524 last successful test
      gapExtend = 2.4, #2.7 last successful test
      factorDiag = 2,
      factorGap = 1,
      binSize = 0.05 # 0.01 last successful test
    ),
    verbose = FALSE
  )
})

saveRDS(fGroups, "fGroups-original.rds")
fGroups <- readRDS("fGroups-original.rds")

# Basic rule based filtering (25353)
fGroups <-
  patRoon::filter(
    fGroups,
    absMinIntensity = 100000, # New - remove noisy peaks
    absMinReplicateAbundance = NULL, # Minimum feature abundance in a replicate group
    relMinReplicateAbundance = 1, # Minimum feature abundance in a replicate group
    relMinReplicates = NULL, # Minimum feature abundance in different replicates
    maxReplicateIntRSD = 0.5, # Maximum relative standard deviation of feature intensities in a replicate group.
    relMinAnalyses = NULL, # Minimum feature abundance in all analyses
    absMinAnalyses = NULL,
    blankThreshold = 3, 
    removeBlanks = FALSE # For validation, maybe don't remove blanks ???
  )

# Save Workspace
rm(fList)
save.image(".RData")

#--------------------------
# NeatMS
#--------------------------

# Export aligned feature groups to .csv for NeatMS analysis
source("https://raw.githubusercontent.com/drewszabo/Rntms/main/create_aligned_table.R")
feature_dataframe <- create_aligened_features(fGroups)

# Run NeatMS analysis (Python/Jupyter)

# Convert NeatMS results to YAML for filtering
source("https://raw.githubusercontent.com/drewszabo/Rntms/main/convert_to_yaml.R")
convert_to_yaml(ntms_results = "neatms_export.csv")

# Filter based on NeatMS prediction model (17314)
fGroups <- patRoon::filter(fGroups,
                           checkFeaturesSession = "model_session.yml",
                           removeBlanks = TRUE) # remove blanks here helped the picking of peaks with mzR here 


# Save Workspace
save.image(".RData")

# -------------------------
# CAMERA Adduct Componentization
# -------------------------

time.camera <- system.time({
  components <- generateComponents(fGroups,
                                   "camera",
                                   ionization = "positive")
})

# Filter based on adduct formation (16027)
fGroups <- selectIons(fGroups, components, "[M+H]+")

# Save Workspace
save.image(".RData")

# -------------------------
# Export NTS Group Results
# -------------------------

resultsfGroups <- patRoon::as.data.table(fGroups, area = TRUE)

# Remove QC for only NTS features
fGroupsNTS <- patRoon::filter(fGroups, rGroups = "QC", negate = TRUE)
fGroupsNTS <- patRoon::filter(fGroups, rAnalysis = "")
resultsfGroupsNTS <- patRoon::as.data.table(fGroupsNTS, area = TRUE, average = TRUE)


# -------------------------
# MS Peak Annotation
# -------------------------

# Set parameters (mz window)
avgFeatParams <- getDefAvgPListParams(clusterMzWindow = 0.005,
                                      topMost = 250
                                      #method = "distance" # default "hclust" uses clustered height
                                      )

precRules <- getDefIsolatePrecParams(maxIsotopes = 4,
                                     mzDefectRange = c(-0.1, 0.1)
                                     )


# Calculate MS and MSMS peak lists from suspect screening

time.mzr <- system.time({
  mslists <- generateMSPeakLists(
    fGroups,
    "mzr",
    maxMSRtWindow = 5,
    precursorMzWindow = 0.2, # +/- 0.2 Da = 0.4 Da
    topMost = NULL,
    avgFeatParams = avgFeatParams,
    avgFGroupParams = avgFeatParams
  )
})


# Filtering only top 95% MSMS peaks based on relative abundance

mslists <- patRoon::filter(
  mslists,
  absMSIntThr = 1000,
  relMSMSIntThr = 0.05, # trying to reduce noise (helped with at least 1)
  withMSMS = TRUE,
  minMSMSPeaks = 1,
  retainPrecursorMSMS = TRUE,
  isolatePrec = precRules, # Issue 87 fixed 24-07-23
)

# save and recover mslists to save read/writes in RData
saveRDS(mslists, "mslists.rds")
rm(mslists)
mslists <- readRDS("mslists.rds")

# Save Workspace
save.image(".RData")

library(patRoon)
library(xcms)
library(IPO)
library(msPurity)

options(patRoon.MP.maxProcs = 4) # Limit max processes in patRoon
register(SerialParam()) # Limit processes to prevent XCMS crashing (RT alignment & grouping)

# -------------------------
# Load analysis information
# -------------------------

# Manually created list from multiple paths
anaInfo <- read.csv("anaInfo.csv", check.names = FALSE)


# -------------------------
# Optimise peak picking and alignment
# -------------------------

# Set datafiles Path from anaInfo
datafiles <- paste0(anaInfo$path, "/", anaInfo$analysis, ".mzML", sep = "")

# Get Default XCMS Parameters
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')

# Set New Optimisation Parameters
peakpickingParameters$min_peakwidth <- c(10, 30) # cannot overlap with max_peakwidth
peakpickingParameters$max_peakwidth <- c(30, 60)
peakpickingParameters$ppm <- c(1, 30)
peakpickingParameters$noise <- 10000

# Set Experimental Parameters
peakpickingParameters$max <- 2

# Run Experiments
time.xcmsSet <- system.time({
  resultPeakpicking <- 
    optimizeXcmsSet(files = datafiles[3:4], # medium level calibrations
                    params = peakpickingParameters, 
                    nSlaves = 1, 
                    subdir = NULL,
                    plot = TRUE)
})

# Show/Save Results

resultPeakpicking$best_settings$result
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset


# Retention Time Optimization

retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$gapExtend <- 2.7
time.RetGroup <- system.time({
  resultRetcorGroup <-
    optimizeRetGroup(xset = optimizedXcmsSetObject, 
                     params = retcorGroupParameters, 
                     nSlaves = 1, 
                     subdir = NULL,
                     plot = TRUE)
})

# Display All Optimisation Settings

writeRScript(resultPeakpicking$best_settings$parameters, 
             resultRetcorGroup$best_settings)

# -------------------------
# Find features from all samples
# -------------------------

# Peak picking (72227)
fList <- findFeatures(anaInfo, "xcms3", 
                      param = xcms::CentWaveParam(peakwidth = c(10, 57),
                                                  ppm = 21.3,
                                                  noise = 10000,
                                                  snthresh = 10,
                                                  mzdiff = 0.01))

saveRDS(fList, "fList-original.rds")
fList <- readRDS("fList-original.rds")

# Group features and align retention time (30958)
fGroups <- groupFeatures(fList, "xcms3", rtalign = TRUE, loadRawData = TRUE,
                         groupParam = xcms::PeakDensityParam(sampleGroups = anaInfo$group,
                                                             minFraction = 0,
                                                             minSamples = 1,
                                                             bw = 0.879999999999999),
                         retAlignParam = xcms::ObiwarpParam(center = 1,
                                                            response = 1,
                                                            gapInit = 0.6112,
                                                            gapExtend = 2.4,
                                                            factorDiag = 2,
                                                            factorGap = 1))

saveRDS(fGroups, "fGroups-original.rds")
fGroups <- readRDS("fGroups-original.rds")

# Basic rule based filtering (26103)
fGroups <-
  patRoon::filter(
    fGroups,
    relMinReplicateAbundance = NULL, # Minimum feature abundance in a replicate group
    relMinReplicates = NULL, # Minimum feature abundance in different replicates
    maxReplicateIntRSD = NULL, # Maximum relative standard deviation of feature intensities in a replicate group.
    relMinAnalyses = NULL, # Minimum feature abundance in all analyses
    blankThreshold = 3,
    removeBlanks = TRUE
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

# Filter based on NeatMS prediction model (5621)
fGroups <- patRoon::filter(fGroups, checkFeaturesSession = "model_session.yml")


# Save Workspace
save.image(".RData")

# -------------------------
# msPurity
# -------------------------

px <- purityX(xset = fGroupsSusp@features@xdata, cores = 1)

test <- fGroupsSusp@features@xdata

# -------------------------
# Suspect Screening
# -------------------------

suspects <- read.csv("Calmix_IE.csv")

suspects[suspects == ""] <- NA

suspects <- suspects %>%
  select(name, SMILES, adduct)

fGroupsSusp <- screenSuspects(fGroups, suspects, mzWindow = 0.010, onlyHits = TRUE)

resultsfGroupsSusp <- patRoon::as.data.table(fGroupsSusp, area = TRUE, average = TRUE)

checkFeatures(fGroupsSusp)

fGroupsSusp <- patRoon::filter(fGroupsSusp, checkFeaturesSession = "checked-features.yml")

# -------------------------
# MS Peak Annotation
# -------------------------

# Set parameters (mz window)
avgFeatParams <- getDefAvgPListParams(clusterMzWindow = 0.002,
                                      topMost = 250,
                                      method = "distance") # default "hclust" uses clustered height

precRules <- getDefIsolatePrecParams(maxIsotopes = 4)


# Calculate MS and MSMS peak lists from suspect screening

time.mzr <- system.time({
  mslists <- generateMSPeakLists(
    fGroupsSusp,
    "mzr",
    maxMSRtWindow = 5,
    precursorMzWindow = 0.2, # +/- 0.2 Da = 0.4 Da
    topMost = NULL,
    avgFeatParams = avgFeatParams,
    avgFGroupParams = avgFeatParams
  )
})


# Filtering only top 99% MSMS peaks based on relative abundance
mslists <- patRoon::filter(mslists,
                           absMSIntThr = 1000,
                           relMSMSIntThr = 0.01,
                           withMSMS = TRUE,
                           minMSMSPeaks = 1,
                           retainPrecursorMSMS = TRUE,
                           isolatePrec = precRules,
                           reAverage = TRUE)

# Save Workspace
save.image(".RData")


# -------------------------
# Formula Generation (SIRIUS)
# -------------------------

time.SIRIUSfor <- system.time({
  formulas <- generateFormulas(
    fGroupsSusp,
    mslists,
    "sirius",
    relMzDev = 5,
    adduct = "[M+H]+",
    elements = siriusElements,
    profile = "orbitrap",
    topMost = 5,
    calculateFeatures = FALSE,
    splitBatches = FALSE,
    projectPath = "log/sirius_formulas/output"
  )
})

# Save Workspace
save.image(".RData")

# -------------------------
# Compound Generation (SIRIUS)
# -------------------------

# Filter for sample group only
fGroupsTest <- patRoon::filter(fGroupsSusp, absMinIntensity = 100000000, mzRange = c(100, 500))

time.SIRIUS <- system.time({
  compoundsSIR <-
    generateCompounds(
      fGroupsSusp,
      mslists,
      "sirius",
      relMzDev = 10,
      adduct = "[M+H]+",
      fingerIDDatabase = "all",
      topMost = 5,
      topMostFormulas = 5,
      profile = "orbitrap",
      elements = siriusElements,
      splitBatches = FALSE,
      cores = 4,
      projectPath = "log/sirius_compounds/output"
    )
})

# Add formula scoring
compoundsSIR <- addFormulaScoring(compoundsSIR, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsSIR <- patRoon::filter(compoundsSIR, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsSIR <- as.data.table(compoundsSIR, fGroups = fGroups)

# Generate toxicity score based on FingerID ???

# Save Workspace
save.image(".RData")


# -------------------------
# Compound Generation (MetFrag)
# -------------------------

time.MetFrag <- system.time({
  compoundsMF <-
    generateCompounds(
      fGroupsSusp,
      mslists,
      "metfrag",
      method = "CL",
      topMost = 5,
      dbRelMzDev = 10,
      fragRelMzDev = 10,
      adduct = "[M+H]+",
      database = "pubchemlite"
    )
})

# Filter for minimum explained peaks and formula score
compoundsMF <- patRoon::filter(compoundsMF, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsMF <- patRoon::as.data.table(compoundsMF, fGroups = fGroups)

# Save Workspace
save.image(".RData")

# -------------------------
# Compound Generation (MassBank)
# -------------------------

mslibrary <- loadMSLibrary("C:/Users/drsz9242/OneDrive - Kruvelab/Drew Szabo/R/MassBank/MassBank_NIST.msp", "msp")

time.MassBank <- system.time({
  compoundsMB <-
    generateCompounds(
      fGroupsSusp,
      mslists,
      "library",
      adduct = "[M+H]+",
      MSLibrary = mslibrary,
      minSim = 0.05,
      absMzDev = 0.010,
      spectrumType = "MS2"
    )
})

# Filter for minimum explained peaks and formula score
compoundsMB <- patRoon::filter(compoundsMB, minExplainedPeaks = 1, topMost = 1)

# Export results as
resultsMB <- patRoon::as.data.table(compoundsMB, fGroups = fGroups)

# Save Workspace
save.image(".RData")


# -------------------------
# reporting (Consensus)
# -------------------------

# Write group summary (total)
resultsExport <- resultsfGroupsSusp %>%
  left_join(resultsMF, by = c("group", "ret"))

report(fGroupsSusp, MSPeakLists = mslists, compounds = compoundsMF)

library(patRoon)
library(xcms)
library(IPO)
library(msPurity)

options(patRoon.MP.maxProcs = 4) # Limit max processes in patRoon
options(future.globals.maxSize = 1000000000)
register(SerialParam()) # Limit processes to prevent XCMS crashing (RT alignment & grouping)

# -------------------------
# Convert raw files to mzML
# -------------------------

convertMSFiles("Raw", "mzML", dirs = TRUE, from = "thermo",
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

# Peak picking ()
fList <- findFeatures(anaInfo, "xcms3",
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
                      verbose = FALSE)

saveRDS(fList, "fList-original.rds")
fList <- readRDS("fList-original.rds")

# Group features and align retention time ()
fGroups <- groupFeatures(fList,
                         "xcms3",
                         rtalign = TRUE,
                         loadRawData = TRUE,
                         groupParam = xcms::PeakDensityParam(sampleGroups = anaInfo$group,
                                                             minFraction = 0,
                                                             minSamples = 1,
                                                             bw = 15, # cranked from 10 due to late eluting big peaks
                                                             binSize = 0.01), # corrected for misaligned m/z in features
                         retAlignParam = xcms::ObiwarpParam(center = 2,
                                                            response = 1,
                                                            gapInit = 0.3, #0.524 last successful test
                                                            gapExtend = 2.4, #2.7 last successful test
                                                            factorDiag = 2,
                                                            factorGap = 1,
                                                            binSize = 0.05), # 0.01 last successful test
                         verbose = FALSE)

saveRDS(fGroups, "fGroups-original.rds")
fGroups <- readRDS("fGroups-original.rds")

# Basic rule based filtering ()
fGroups <-
  patRoon::filter(
    fGroups,
    absMinReplicateAbundance = NULL, # Minimum feature abundance in a replicate group
    relMinReplicateAbundance = NULL, # Minimum feature abundance in a replicate group
    relMinReplicates = NULL, # Minimum feature abundance in different replicates
    maxReplicateIntRSD = NULL, # Maximum relative standard deviation of feature intensities in a replicate group.
    relMinAnalyses = NULL, # Minimum feature abundance in all analyses
    absMinAnalyses = 2,
    blankThreshold = NULL, # For validation, maybe don't remove blanks ???
    removeBlanks = FALSE
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
fGroups <- patRoon::filter(fGroups,
                           checkFeaturesSession = "model_session.yml",
                           removeBlanks = TRUE) # remove blanks here helped the picking of peaks with mzR here 


# Save Workspace
save.image(".RData")

# -------------------------
# msPurity
# -------------------------

px <- purityX(xset = fGroupsSusp@features@xdata, cores = 1)

test <- fGroupsSusp@features@xdata






components <- generateComponents(
  fGroups,
  "camera",
  ionization = "positive")

fGroups <- selectIons()

# -------------------------
# Suspect Screening
# -------------------------

suspects <- read.csv("Calmix_IE.csv")

suspects[suspects == ""] <- NA

suspects <- suspects %>%
  select(name, SMILES, adduct) %>% # removed rt variable to include all features with relevant mz
  drop_na(name)

fGroupsSusp <- screenSuspects(fGroups, suspects, mzWindow = 0.025, onlyHits = TRUE)

checkFeatures(fGroupsSusp)

fGroupsSusp <- patRoon::filter(fGroupsSusp, checkFeaturesSession = "checked-features.yml")

resultsfGroupsSusp <- patRoon::as.data.table(fGroupsSusp, area = TRUE, average = TRUE)


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


# Filtering only top 99% MSMS peaks based on relative abundance

mslists <- patRoon::filter(mslists,
                           absMSIntThr = 1000,
                           relMSMSIntThr = 0.05, # trying to reduce noise (helped with at least 1)
                           withMSMS = TRUE,
                           minMSMSPeaks = 1,
                           retainPrecursorMSMS = TRUE,
                           isolatePrec = precRules, # Issue 87 fixed 24-07-23
                           )

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
    elements = "CHONPSFClBr",
    profile = "orbitrap",
    topMost = 10,
    calculateFeatures = FALSE,
    splitBatches = FALSE
  )
})

# Save Workspace
save.image(".RData")



# -------------------------
# Compound Generation (SIRIUS)
# -------------------------

time.SIRIUS <- system.time({
  compoundsSIR <-
    generateCompounds(
      fGroupsSusp,
      mslists,
      "sirius",
      relMzDev = 5,
      adduct = "[M+H]+",
      formulaDatabase = "pubchem",
      topMost = 5,
      topMostFormulas = 10, # from 5 - hopefully increase the number of form used to calculate structures
      profile = "orbitrap",
      splitBatches = FALSE,
      cores = 4,
      elements = "CHONPSFClBr",
      extraOptsFormula = "--ppm-max-ms2=50",
      verbose = TRUE
    )
})



# Add formula scoring
compoundsSIR <- addFormulaScoring(compoundsSIR, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsSIR <- patRoon::filter(compoundsSIR, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsSIR <- patRoon::as.data.table(compoundsSIR, fGroups = fGroups)

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
      dbRelMzDev = 5,
      fragAbsMzDev = 0.02, # changed from 5 ppm (relative) to equal MassBank
      adduct = "[M+H]+",
      database = "pubchemlite",
      maxCandidatesToStop = 2500 # resource intensive - consider using PubChemLite to reduce #candidates
    )
})

# Add formula scoring
compoundsMF <- addFormulaScoring(compoundsMF, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsMF <- patRoon::filter(compoundsMF, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsMF <- patRoon::as.data.table(compoundsMF, fGroups = fGroups)

# Save Workspace
save.image(".RData")

# -------------------------
# Compound Generation (MassBank)
# -------------------------

mslibrary <- loadMSLibrary("C:/Data/MassBank/MassBank_NIST.msp", "msp")

simParam <- getDefSpecSimParams(
  absMzDev = 0.02 # 20 mDa difference for MS2 spectra
  ) # https://rickhelmus.github.io/patRoon/reference/specSimParams.html

time.MassBank <- system.time({
  compoundsMB <-
    generateCompounds(
      fGroupsSusp,
      mslists,
      "library",
      adduct = "[M+H]+",
      MSLibrary = mslibrary,
      minSim = 0.50,
      absMzDev = 0.05,
      spectrumType = "MS2",
      checkIons = "adduct",
      specSimParams = simParam # increase bin size
    )
})

# Add formula scoring
compoundsMB <- addFormulaScoring(compoundsMB, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsMB <- patRoon::filter(compoundsMB, minExplainedPeaks = 1, topMost = 1)

# Export results as
resultsMB <- patRoon::as.data.table(compoundsMB, fGroups = fGroups)

# Save Workspace
save.image(".RData")


# -------------------------
# MS2Tox - Toxicity Predictions
# -------------------------

library(MS2Tox)

logPath = "log/sirius_compounds/sirius-batch_1-[M+H]+.txt"
logFile <- readr::read_file(logPath)

folderwithSIRIUSfiles <- stringr::word(logFile, 7,7) # dumb location setting, may change with version
folderwithSIRIUSfiles <- gsub("\\", "/", folderwithSIRIUSfiles, fixed = TRUE)

UnZip_SIRIUS5(folderwithSIRIUSfiles)

resultsMS2Tox  <- FishLC50Prediction(folderwithSIRIUSfiles, "static")

# -------------------------
# reporting (Consensus)
# -------------------------

# Write group summary (total)
resultsExport <- resultsfGroupsSusp %>%
  left_join(resultsMF, by = c("group", "ret")) %>%
  left_join(resultsSIR, by = c("group", "ret")) %>%
  left_join(resultsMB, by = c("group", "ret"))

write.csv(resultsExport, "resultsID.csv")

compounds <- consensus(compoundsMF, compoundsSIR)

# SIRIUS Report
report(fGroupsSusp, MSPeakLists = mslists, formulas = formulas, compounds = compoundsSIR)

# MetFrag Report
report(fGroupsSusp, MSPeakLists = mslists, formulas = formulas, compounds = compoundsMF)

# MassBank Report
report(fGroupsSusp, MSPeakLists = mslists, formulas = formulas, compounds = compoundsMB)


plotChroms(fGroups[, "M297_R498_10701"], # only plot all features of first group
           colourBy = "rGroup") # and mark them individually per replicate group


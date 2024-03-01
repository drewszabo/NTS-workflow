library(patRoon)

# -------------------------
# Formula Generation (SIRIUS)
# -------------------------

time.SIRIUSfor <- system.time({
  formulas <- generateFormulas(
    fGroups,
    MSPeakLists = readRDS("mslists.rds"),
    "sirius",
    relMzDev = 5,
    adduct = "[M+H]+",
    database = "pubchem",
    topMost = 10,
    profile = "orbitrap",
    splitBatches = FALSE,
    cores = 4,
    elements = "CHONPSFClBr",
    extraOptsFormula = "--ppm-max-ms2=50",
    calculateFeatures = FALSE,
    verbose = TRUE
  )
})

# Filter formulas
formulas <- patRoon::filter(formulas, minExplainedPeaks = 2)

# Save Workspace
save.image(".RData")



# -------------------------
# Compound Generation (SIRIUS)
# -------------------------

time.SIRIUS <- system.time({
  compoundsSIR <-
    generateCompounds(
      fGroups,
      MSPeakLists = readRDS("mslists.rds"),
      "sirius",
      relMzDev = 5,
      adduct = "[M+H]+",
      formulaDatabase = "pubchem",
      topMost = 5,
      topMostFormulas = 10, # from 5 - hopefully increase the number of form used to calculate structures
      profile = "orbitrap",
      splitBatches = FALSE,
      cores = 4, # important for efficient processing
      elements = "CHONPSFClBr",
      extraOptsFormula = "--ppm-max-ms2=50",
      verbose = TRUE
    )
})

# Get folder path from patRoon log
folderwithSIRIUSfiles <- "log/sirius_compounds/sirius-batch_1-[M+H]+.txt"
folderwithSIRIUSfiles <- readr::read_file(folderwithSIRIUSfiles)
folderwithSIRIUSfiles <- stringr::word(folderwithSIRIUSfiles, 6,6) # dumb location setting, may change with version
folderwithSIRIUSfiles <- gsub("\\", "/", folderwithSIRIUSfiles, fixed = TRUE)

# Copy files to home directory
file.copy(folderwithSIRIUSfiles, "log/sirius_compounds",
          recursive = TRUE)
folderwithSIRIUSfiles <- dir("log/sirius_compounds", full.names = TRUE)[2]

# Backup
saveRDS(compoundsSIR, "compoundsSIR.rds")

# Add formula scoring
compoundsSIR <- addFormulaScoring(compoundsSIR, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsSIR <- patRoon::filter(compoundsSIR, minExplainedPeaks = 2, topMost = 1) # consider filtering for RTI instead of Top 1

# Export results as
resultsSIR <- patRoon::as.data.table(compoundsSIR, fGroups = fGroups) %>%
  select(group, ret, compoundName, score, SMILES, identifier, neutral_formula) %>%
  rename(compoundNameSIR = compoundName,
         scoreSIR = score,
         SMILESSIR = SMILES,
         cidSIR = identifier,
         formulaSIR = neutral_formula) %>%
  mutate(cidSIR = str_split_i(cidSIR, ";", 1),
         ret = ret/60)

# Save Workspace
save.image(".RData")


# -------------------------
# Compound Generation (MetFrag)
# -------------------------

time.MetFrag <- system.time({
  compoundsMF <-
    generateCompounds(
      fGroups,
      MSPeakLists = readRDS("mslists.rds"),
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

# Backup
saveRDS(compoundsMF, "compoundsMF.rds")

# Add formula scoring
compoundsMF <- addFormulaScoring(compoundsMF, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsMF <- patRoon::filter(compoundsMF, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsMF <- patRoon::as.data.table(compoundsMF, fGroups = fGroups)

# Isolate MoNA scores
resultsMoNA <- resultsMF %>%
  dplyr::filter(individualMoNAScore > 0.5) %>%
  dplyr::select(group, ret, compoundName, individualMoNAScore, SMILES, InChIKey, neutral_formula) %>%
  dplyr::rename(compoundNameMoNA = compoundName,
         scoreMoNA = individualMoNAScore,
         SMILESMoNA = SMILES,
         InChIKeyMoNA = InChIKey,
         formulaMoNA = neutral_formula)

resultsMF <- resultsMF %>%
  dplyr::select(group, ret, compoundName, score, SMILES, InChIKey, identifier, neutral_formula) %>%
  dplyr::rename(compoundNameMF = compoundName,
         scoreMF = score,
         SMILESMF = SMILES,
         InChIKeyMF = InChIKey,
         cidMF = identifier,
         formulaMF = neutral_formula) %>%
  dplyr::mutate(cidMF = str_split_i(cidMF, " ", 1),
                ret = ret/60)



# Save Workspace
save.image(".RData")

# -------------------------
# Compound Generation (MassBank)
# -------------------------

simParam <- getDefSpecSimParams(
  absMzDev = 0.02 # 20 mDa difference for MS2 spectra
) # https://rickhelmus.github.io/patRoon/reference/specSimParams.html

time.MassBank <- system.time({
  compoundsMB <-
    generateCompounds(
      fGroups,
      MSPeakLists = readRDS("mslists.rds"),
      "library",
      adduct = "[M+H]+",
      MSLibrary = loadMSLibrary("C:/Data/MassBank/MassBank_NIST.msp", "msp"),
      minSim = 0.50,
      absMzDev = 0.05,
      spectrumType = "MS2",
      checkIons = "adduct",
      specSimParams = simParam # increase bin size
    )
})

# Backup
saveRDS(compoundsMB, "compoundsMB.rds")

# Add formula scoring
compoundsMB <- addFormulaScoring(compoundsMB, formulas, updateScore = TRUE)

# Filter for minimum explained peaks and formula score
compoundsMB <- patRoon::filter(compoundsMB, minExplainedPeaks = 2, topMost = 1)

# Export results as
resultsMB <- patRoon::as.data.table(compoundsMB, fGroups = fGroups) %>%
  select(group, compoundName, libMatch, SMILES, InChIKey, neutral_formula) %>%
  rename(compoundNameMB = compoundName,
         scoreMB = libMatch,
         SMILESMB = SMILES,
         InChIKeyMB = InChIKey,
         formulaMB = neutral_formula)

# Save Workspace
save.image(".RData")



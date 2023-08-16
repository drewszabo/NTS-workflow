# NTS-workflow

## Description
Utilizing patRoon and other tools to perform non-target analysis of high-resolution mass spectrometry data

## Scope

This workflow was optimized for LC-ESI-HRMS data with data-dependant acquisition (DDA). Separation was achieved with a 25 minute elution (A = 0.1% fomic acid, B = acetonitrile) and a reverse phase C18 column (100 mm x 2.1 mm). Mass spectrum were acquired with Thermo Q Exactive (Thermo) in ESI+ mode.

## Dependancies
```
library(tidyverse)
library(patRoon)
library(xcms) # installed by patRoon but often helps if loaded seperately
library(IPO)
```

## Workflow

### 1. Feature Extraction

#### IPO

The ‘Isotopologue Parameter Optimization’ (`IPO`) tool was used to optimise the `xcms` peak picking and alignment parameters.

Not used in final process. IPO was not able to correctly optimise peak picking or grouping for my dataset. Perhaps it will be useful under different LC-HRMS conditions.

<details>
  <summary>Show code</summary>

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("IPO")

# Get Default XCMS Parameters
peakpickingParameters <- getDefaultXcmsSetStartingParams('centWave')

# Set New Optimisation Parameters
peakpickingParameters$min_peakwidth <- c(6, 18)
peakpickingParameters$max_peakwidth <- c(30, 90)
peakpickingParameters$ppm <- c(5,40)
peakpickingParameters$mzdiff <- c(-0.01, -0.001)
peakpickingParameters$snthresh <- c(3, 17)
peakpickingParameters$noise <- c(0, 5000)

# Run Experiments
time.xcmsSet <- system.time({ # measuring time
  resultPeakpicking <- 
    optimizeXcmsSet(files = datafiles[1:6], 
                    params = peakpickingParameters, 
                    nSlaves = 1, 
                    subdir = NULL,
                    plot = TRUE)
})

# Show/Save Results
resultPeakpicking$best_settings$result
optimizedXcmsSetObject <- resultPeakpicking$best_settings$xset

# Retention Time / Alignment Optimisation
retcorGroupParameters <- getDefaultRetGroupStartingParams()
retcorGroupParameters$profStep <- 1
retcorGroupParameters$gapExtend <- 2.7

time.RetGroup <- system.time({ # measuring time
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
```

</details>

Libiseller, G., Dvorzak, M., Kleb, U., Gander, E., Eisenberg, T., Madeo, F., Neumann, S., Trausinger, G., Sinner, F., Pieber, T. & Magnes, C. 2015. IPO: a tool for automated optimization of XCMS parameters. BMC Bioinformatics, 16, 118. https://doi.org/10.1186/s12859-015-0562-8

#### XCMS (via patRoon)

Using optimised parameters from IPO, all features were extracted from the `.mzML` files. Note: current version implementation of XCMS is not always compatable with multithreded systems. Single core analysis may be required: `register(SerialParam())`

<details>
  <summary>Show code</summary>

```
# Extract all features
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

# Perform feature alignment
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
```

</details>

Smith, C.A., Want, E.J., O'Maille, G., Abagyan, R. & Siuzdak, G. 2006. XCMS:  Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. Analytical Chemistry, 78, 779-787. https://doi.org/10.1021/ac051437y

#### Filtering

Feature groups are filtered for replicate and blank intensity/abundance.

<details>
  <summary>Show code</summary>

```
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
```

</details>

#### NeatMS

A table of feature parameters is generated and exported for implementation in Python (see Documentation from authors). See [drewszabo/neatms.export](https://www.github.com/drewszabo/ntms.export) for updated code to produce data table that is compatable with this workflow.

<details>
  <summary>Show code</summary>

```
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

```

Gloaguen, Y., Kirwan, J.A. & Beule, D. 2022. Deep Learning-Assisted Peak Curation for Large-Scale LC-MS Metabolomics. Analytical Chemistry, 94, 4930-4937. https://doi.org/10.1021/acs.analchem.1c02220

#### msPurity

The quality of MS2 spectra can be evaluated by assessing the purity of the MS1 peaks present in the isolation window. If chimeric peaks are detected, the score will be reduced and the feature can be ommited from further analysis, due to poor MS2 spectrum quality. The authors recommend a minimum score of 0.5 to continue with peak annotation and identification.
  
  It is hypothesised that this may also reduce the number of "noisy" EIC, as the ratio of MS1 peaks could be reduced if there is low abundance of the selected peak.

<details>
  <summary>Show code</summary>

Code not yet implemented or tested. -DS

</details>

Lawson, T.N., Weber, R.J.M., Jones, M.R., Chetwynd, A.J., Rodrı́guez-Blanco, G., Di Guida, R., Viant, M.R. & Dunn, W.B. 2017. msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry-Based Fragmentation in Metabolomics. Analytical Chemistry, 89, 2432-2439. https://doi.org/10.1021/acs.analchem.6b04358

### 2. MS Peak Annotation
  
  #### mzR (via patRoon)
  
  Default mzR parameters to calculate the average peak lists were changed to better suit current workflows. May require further optimisation - eg. topMost = 250 significantly increases compute time.
  
  Filtering includes precursor isolation and MS2 abundance to clean spectra and significantly reduce object size.
  
  <details>
  <summary>Show code</summary>

```
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
```

</details>
    
### 3. Compound Annotation
    
For each of the following library and in-silico matching, a minimum of 2 explained peaks is used as a filter to increase the confidence of annotation. In general, the score can be artificially inflated for a compound that only matches with one (or 0) fragment, due to limitations of the cosine (dot-product) scoring technique.
    
#### MassBank (via patRoon)
    
Requires latest database `.msp` download from MassBank repo
    
<details>
  <summary>Show code</summary>

```
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

# Filter for minimum explained peaks and formula score
compoundsMB <- patRoon::filter(compoundsMB, topMost = 1, minExplainedPeaks = 2)

# Export results as
resultsMB <- patRoon::as.data.table(compoundsMB, fGroups = fGroups)
```

</details>

#### SIRIUS CSIFingerID
  
Currently working with SIRIUS v5.6.3. Limiting the cores may not be necessary but can help if you need to use your computer while processing the data. Using the projectPath variable is also not necessary, although I find it useful to browse the raw SIRIUS results for troubleshooting.
  
<details>
<summary>Show code</summary>

```
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
                  
  # Filter for minimum explained peaks and SIRIUS score
  compoundsSIR <- patRoon::filter(compoundsSIR, topMost = 1, minExplainedPeaks = 2)

  # Export results as
  resultsSIR <- patRoon::as.data.table(compoundsSIR, fGroups = fGroups)
                  
  ```

</details>
  
#### MetFrag
  
  Current issues with mass accuracy have required relatively high ppm deviations.
  
  <details>
  <summary>Show code</summary>

```
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

# Filter for minimum explained peaks and formula score
compoundsMF <- patRoon::filter(compoundsMF, topMost = 1, minExplainedPeaks = 2)

# Export results as
resultsMF <- patRoon::as.data.table(compoundsMF, fGroups = fGroups)
```

</details>

# NTS-workflow

## Description
Utilizing patRoon and other tools to perform non-target analysis of high-resolution mass spectrometry data

## Scope

This workflow was optimized for LC-ESI-HRMS data with data-dependant acquisition (DDA). Separation was achieved with a 25 minute elution (A = 0.1% fomic acid, B = acetonitrile) and a reverse phase C18 column (100 mm x 2.1 mm). Mass spectrum were acquired with Thermo Q Exactive (Thermo) in ESI+ mode.

## Workflow

### 1. Feature Extraction

#### IPO

The ‘Isotopologue Parameter Optimization’ (`IPO`) tool was used to optimise the `xcms` peak picking and alignment parameters.

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
fList <- patRoon::findFeatures(
  anaInfo,
  "xcms3",
  param = xcms::CentWaveParam(
    ppm = 18.5,
    mzdiff = -0.0145,
    prefilter = c(3, 100),
    snthresh = 4.4,
    peakwidth = c(9, 76),
    noise = 7500
  )
)

# Perform feature alignment
fGroups <-
  groupFeatures(
    fList,
    "xcms3",
    rtalign = TRUE,
    loadRawData = TRUE,
    groupParam = xcms::PeakDensityParam(
      sampleGroups = anaInfo$group,
      minFraction = 0,
      minSamples = 1,
      bw = 0.87999
    ),
    retAlignParam = xcms::ObiwarpParam(
      gapInit = 0.8416,
      gapExtend = 2.7,
      factorDiag = 2,
      factorGap = 1,
      response = 1,
      centerSample = 3
    )
  )
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
    relMinReplicateAbundance = 1,
    maxReplicateIntRSD = 0.50,
    blankThreshold = 3,
    removeBlanks = TRUE
  )
```

</details>

#### NeatMS

A table of feature parameters is generated and exported for implementation in Python (see Documentation from authors). See [drewszabo/neatms.export](https://www.github.com/drewszabo/ntms.export) for updated code to produce data table that is compatable with this workflow.

Gloaguen, Y., Kirwan, J.A. & Beule, D. 2022. Deep Learning-Assisted Peak Curation for Large-Scale LC-MS Metabolomics. Analytical Chemistry, 94, 4930-4937. https://doi.org/10.1021/acs.analchem.1c02220

### 2. MS Peak Annotation
  
  #### mzR (via patRoon)
  
  Default mzR parameters to calculate the average peak lists were changed to better suit current workflows. May require further optimisation - eg. topMost = 250 significantly increases compute time.
  
  Filtering includes precursor isolation and MS2 abundance to clean spectra and significantly reduce object size.
  
  <details>
  <summary>Show code</summary>

```
# Set parameters (mz window)
avgFeatParams <- getDefAvgPListParams(
  clusterMzWindow = 0.002,
  topMost = 250,
  minIntensityPre = 500,
  minIntensityPost = 1000,
  method = "hclust",
  pruneMissingPrecursorMS = TRUE,
  retainPrecursorMSMS = TRUE
)


# Calculate MS and MSMS peak lists from suspect screening
mslists <- generateMSPeakLists(
  fGroups,
  "mzr",
  maxMSRtWindow = 15,
  precursorMzWindow = 0.4,
  topMost = NULL,
  avgFeatParams = avgFeatParams,
  avgFGroupParams = avgFeatParams
)


# Filtering only top 99% MSMS peaks based on relative abundance
mslists <- patRoon::filter(
  mslists,
  absMSIntThr = 1000,
  relMSMSIntThr = 0.01,
  withMSMS = TRUE,
  minMSMSPeaks = 1,
  retainPrecursorMSMS = TRUE,
  isolatePrec = TRUE
)
```

</details>
    
### 3. Compound Annotation
    
#### MassBank (via patRoon)
    
Requires latest database `.msp` download from MassBank repo
    
<details>
  <summary>Show code</summary>

```
compoundsMB <-
  generateCompounds(
    fGroupsSusp,
    mslists,
    "library",
    adduct = "[M+H]+",
    MSLibrary = mslibrary,
    minSim = 0.05,
    absMzDev = 0.01,
    spectrumType = "MS2"
  )

# Filter for minimum explained peaks and formula score
compoundsMB <- patRoon::filter(compoundsMB, topMost = 1)

# Export results as
resultsMB <- patRoon::as.data.table(compoundsMB, fGroups = fGroups)
```

</details>
  
#### MetFrag
  
  Current issues with mass accuracy have required relatively high ppm deviations.
  
  <details>
  <summary>Show code</summary>

```
compoundsMF <-
  generateCompounds(
    fGroupsSusp,
    mslists,
    "metfrag",
    method = "CL",
    topMost = 5,
    dbRelMzDev = 25,
    fragRelMzDev = 25,
    adduct = "[M+H]+",
    database = "pubchemlite"
  )

# Filter for minimum explained peaks and formula score
compoundsMF <- patRoon::filter(compoundsMF, topMost = 1)

# Export results as
resultsMF <- patRoon::as.data.table(compoundsMF, fGroups = fGroups)
```

</details>

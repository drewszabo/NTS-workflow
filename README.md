# NTS-workflow

## Description
Utilizing patRoon and other tools to perform non-target analysis of high-resolution mass spectrometry data

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

Smith, C.A., Want, E.J., O'Maille, G., Abagyan, R. & Siuzdak, G. 2006. XCMS:  Processing Mass Spectrometry Data for Metabolite Profiling Using Nonlinear Peak Alignment, Matching, and Identification. Analytical Chemistry, 78, 779-787. 10.1021/ac051437y

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

A table of feature parameters is exported and 

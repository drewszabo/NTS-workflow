export_neatms <- function(fGroups = fGroupsXCMS, csv_path = "aligned_feature_table.csv"){

# From NeatMS Documentation
# Feature table export from XCMS (aligned peaks with gapfilling)


# The first part of the code is the same as for unaligned peaks, you can jump to the feature information addition 

# This code assumes that the xdata variable corresponds 
# to the XCMSnExp object that contains the detected peaks 

# Load dplyr (required for left_join())
library(tidyverse)
# Load tibble (required for rownames_to_column())
library(tibble)
# Load XCMS
library(xcms)

#--------------
# Extract XCMS object from fGroups (Drew)
#--------------
xdata_filled <- fGroups@xdata

# Create dataframe containing adjusted retention times 
df_chrom_peaks_fill <- as.data.frame(chromPeaks(xdata_filled))
# Create dataframe with raw retention times 
feature_dataframe <- as.data.frame(chromPeaks(dropAdjustedRtime(xdata_filled)))

# Get the peaks that have been recovered by the gapfilling step
df_filled_peaks <- df_chrom_peaks_fill[!row.names(df_chrom_peaks_fill) %in% c(rownames(feature_dataframe)),]

# Add them to the raw retention time dataframe (filled peaks do not have raw retention times so we use the adjusted ones)
feature_dataframe <- dplyr::bind_rows(feature_dataframe, df_filled_peaks)

# Rename the retention time columns of the adjusted rt dataframe
df_chrom_peaks_fill <- df_chrom_peaks_fill %>%
  dplyr::rename(rt_adjusted = rt, rtmin_adjusted = rtmin, rtmax_adjusted = rtmax)

# Add the adjusted rt columns to the dataframe containing raw rt
# To have a dataframe containing both raw and adjusted rt information
feature_dataframe <- dplyr::left_join(rownames_to_column(feature_dataframe), rownames_to_column(df_chrom_peaks_fill[,c("rt_adjusted","rtmin_adjusted","rtmax_adjusted")]), by="rowname")

# Remove the rownames as we won't need them
feature_dataframe$rowname <- NULL

# Retrieve the sample names and store them as a dataframe
sample_names_df <- as.data.frame(MSnbase::sampleNames(xdata_filled))

#------------------
# Add sample names manually from patRoon data (Drew)
#------------------
sample_names_df <- as.data.frame(fGroups@analysisInfo$analysis)
sample_names_df$`fGroups@analysisInfo$analysis` <- paste0(sample_names_df$`fGroups@analysisInfo$analysis`, ".mzML")

# Rename the unique column "sample_name"
colnames(sample_names_df) <- c("sample_name")

# Generate the correct sample ids for matching purposes
# MSnbase sampleNames() function returns sample names ordered by their ids
sample_names_df$sample <- seq.int(nrow(sample_names_df))

# Attach the sample names to the main dataframe by matching ids (sample column)
feature_dataframe <- dplyr::left_join(feature_dataframe,sample_names_df, by="sample")

### Feature information addition ###

# Here we will bring the feature alignment information stored in the XCMSnExp object to the dataframe that we have already created

featuresDef <- featureDefinitions(xdata_filled)
featuresDef_df = data.frame(featuresDef)

#--------------------
# Added as.numeric() to fix the column_index assignment (Drew)
#--------------------

# Adjust variable
# Only keep the information we need (column named 'peakidx')
# Get the index of the peakidx column
column_index <- as.numeric(which(colnames(featuresDef_df)=="peakidx"))
# Extract the peakidx column
features_df <- as.list(featuresDef_df$peakidx)

for(x in 1:length(features_df)){
  length(features_df[[x]]) <- length(fGroups@groups[[x]])
}

features_df = data.frame(features_df)

features_df = data.frame(t(features_df))

# Rename the column
peak_colummn_name <- colnames(features_df)
features_df = dplyr::rename(features_df, "peak_id"=all_of(peak_colummn_name))

#------------------
# Add feature names from patRoon to features_df (Drew)
#------------------

feature_names <- colnames(fGroups@groups)

features_df <- cbind(feature_id = feature_names, features_df)

features_df <- features_df %>%
  pivot_longer(!feature_id, values_to = "peak_id") %>%
  select(feature_id, peak_id) %>%
  drop_na(peak_id)


# We'll use data.table for the next step
require(data.table)

# Get all the peak_id for each feature_id
features_df <- data.table(features_df)
features_df = features_df[, list(peak_id = unlist(peak_id)), by=feature_id]

# Bring the feature_id to the original peak dataframe
feature_dataframe = cbind(peak_id= row.names(feature_dataframe),feature_dataframe)
feature_dataframe$peak_id = as.character(feature_dataframe$peak_id)
features_df$peak_id = as.character(features_df$peak_id)
feature_dataframe = left_join(feature_dataframe, features_df, by="peak_id")

# Note: The dataframe contains an extra column called peak_id, but this won't affect NeatMS and will simply be ignored (as would any other column not present in the list above).

# Export the data as csv. 
# Note: Set row.names to FALSE as NeatMS does not need them
write.csv(feature_dataframe, csv_path, row.names = FALSE)

return(feature_dataframe)

}




import_neatms <- function(ntms_results = "ntms_export.csv", anaInfo = anaInfo, feature_dataframe = feature_dataframe, yaml_path = "model_session.yml"){

# Import results
neatms_export <- read.csv(ntms_results)

# Select relevant columns from feature_dataframe
feature_dataframe <- feature_dataframe %>%
  dplyr::select(feature_id, into)

    if (missing(feature_dataframe)) stop("feature_dataframe must be provided")
  if (missing(anaInfo)) stop("anaInfo must be provided")

# Rename and transform columns in NeatMS export
neatms_export <- neatms_export %>%
  dplyr::rename(
    feature_id = feature.id,
    mz = m.z,
    maxo = height,
    into = area,
    rt = retention.time,
    analysis = sample
  ) %>%
  dplyr::mutate(rt = rt * 60)  # Convert retention time from minutes to seconds

# Join metadata from anaInfo and feature_dataframe
neatms_export <- neatms_export %>%
  dplyr::left_join(dplyr::select(anaInfo, group, analysis), by = "analysis") %>%
  dplyr::left_join(feature_dataframe, by = "into")

write.csv(neatms_export, "neatms_export_aligned.csv", row.names = FALSE)
  
# Summarize feature quality by group
neatms_export <- neatms_export %>%
  dplyr::group_by(feature_ID, group) %>%
  dplyr::summarise(
    quality = any(label == "High_quality"),
    feature = dplyr::first(feature_id),
    .groups = "drop"
  )



# Filter features with <1 "high-quality"
removeFully <- neatms_export %>%
  dplyr::group_by(feature) %>%
  dplyr::filter(sum(quality)/length(quality) == 0) %>%
  dplyr::distinct(feature) %>%
  dplyr::select(feature) %>%
  dplyr::rename(removeFully = feature)


# Create YAML file for removal of noisy features
library(yaml)

list <- as.list(removeFully)

# Set the version field
version <- 1.0 # Should be written without quotation marks

# Set the type field
type <- "featureGroups"

# Set the removePartially field
removePartially <- list()
featureGroups <- list()

# Add the version, type, and removePartially fields to the list
list$removePartially <- removePartially
list$version <- version
list$type <- type

featureGroups <- patRoon::as.data.table(fGroups) %>%
  dplyr::select(group, ret, mz) %>%
  dplyr::filter(group %in% removeFully$removeFully)
featureGroups <- list(featureGroups = split(replace(featureGroups, "group", NULL), featureGroups$group))
list <- append(list, featureGroups)

# Convert the list to a YAML object
yaml_object <- as.yaml(list)

# Write the YAML object to a file
write(yaml_object, file = yaml_path)

}





generateGNPSInfo <- function(fGroups, mslists, path) {
  # Generate feature table
  featuresDef <- xcms::featureDefinitions(fGroups@xdata)
  featuresIntensities <- xcms::featureValues(fGroups@xdata, value = "into")
  peakID <- as.character(colnames(fGroups@groups))
  dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
  dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]
  dataTable$peak_ID <- peakID
  
  write.table(dataTable, paste0(path, "output_gnps.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  resultsmslists <- patRoon::as.data.table(mslists)
  groups <- unique(resultsmslists$group)
  
  fileConn <- file(paste0(path, "output_gnps.mgf"), "w")
  writeLines(paste0("COM=Exported by ", Sys.getenv("USERNAME"), " on ", Sys.time()), fileConn)
  
  for (peak_ID in groups) {
    feature_list <- dataTable[dataTable$peak_ID == peak_ID, ]
    peaks <- resultsmslists[group == peak_ID]
    MS1peaks <- peaks[type == "MS", .(mz, intensity)]
    MS2peaks <- peaks[type == "MSMS", .(mz, intensity)]
    
    ret <- round(feature_list$rtmed, 2)
    mz <- round(feature_list$mzmed, 4)
    numPeaksMS1 <- nrow(MS1peaks)
    numPeaksMS2 <- nrow(MS2peaks)
    
    # MS1 Block
    writeLines(c("BEGIN IONS",
                 paste0("FEATURE_ID=", feature_list$Row.names),
                 paste0("PEAK_ID=", peak_ID),
                 "MSLEVEL=1",
                 paste0("RTINSECONDS=", ret),
                 paste0("PEPMASS=", mz),
                 paste0("SCANS=", as.integer(sub("^FT", "", feature_list$Row.names))),
                 paste0("Num peaks=", numPeaksMS1),
                 apply(MS1peaks, 1, function(x) paste(round(x[1], 6), round(x[2], 0), sep = " ")),
                 "END IONS", ""), fileConn)
    
    # MS2 Block
    writeLines(c("BEGIN IONS",
                 paste0("FEATURE_ID=", feature_list$Row.names),
                 paste0("PEAK_ID=", peak_ID),
                 "MSLEVEL=2",
                 paste0("RTINSECONDS=", ret),
                 paste0("PEPMASS=", mz),
                 paste0("SCANS=", as.integer(sub("^FT", "", feature_list$Row.names))),
                 paste0("Num peaks=", numPeaksMS2),
                 apply(MS2peaks, 1, function(x) paste(round(x[1], 6), round(x[2], 0), sep = " ")),
                 "END IONS", ""), fileConn)
  }
  
  close(fileConn)
}



export_pfascreen <- function(mslists = mslists, path = "results"){

  # Export feature list
ms1int <- patRoon::as.data.table(mslists) %>%
  dplyr::filter(type == "MS") %>%
  group_by(group) %>%
  slice_head(n = 3) %>%
  mutate(ID = row_number()) %>%
  ungroup() %>%
  dplyr::select(group, ID, intensity) %>%
  pivot_wider(names_from = ID, values_from = intensity)

ms1mz <- resultsfGroups %>%
  dplyr::select(group, mz, ret) %>%
  dplyr::rename(rt = ret) %>%
  dplyr::mutate(rt = rt/60) # rt in min, not seconds

pfascreen <- ms1mz %>%
  left_join(ms1int, by = "group") %>%
  drop_na("1") %>%
  mutate(index = row_number() - 1)

colnames(pfascreen) <- c("group", "mz",	"rt",	"mz_area",	"mz+1_area",	"mz+2_area", "index")

return(pfascreen)

write.csv(pfascreen, paste0(path, "/pfascreen.csv", sep = ""), row.names = FALSE, na = "")
  
}





import_pfascreen <- function(folder = "results/PFAScreen-main"){

# List all files in the working directory
file <- list.files(path = folder, pattern = "Results_PFAScreen.*\\.csv", recursive = TRUE, full.names = TRUE)

resultsPFAScreen <- read.csv(file) %>%
  dplyr::rename(index = X) %>%
  dplyr::select(-(2:14))

resultsPFAScreen <- pfascreen %>%
  dplyr::select(group, index) %>%
  left_join(resultsPFAScreen, by = "index") %>%
  dplyr::filter(unique_homologues > 1) %>%
  dplyr::select(-index)

  return(resultsPFAScreen)

}               
               

padelDesc <- function(df, jar_path = "descriptor-cli-0.1a-SNAPSHOT-all.jar") {
  stopifnot("SMILES" %in% names(df))
  base_command <- paste("java -jar", jar_path)
  
  calc_by_smiles <- function(smiles) {
    command_final <- paste(base_command, shQuote(smiles))
    
    java_output <- tryCatch(
      system(command_final, intern = TRUE),
      error = function(e) return(character())
    )
    
    if (length(java_output) == 0 || any(grepl("Exception|Error|Failed", java_output))) {
      message("Descriptor generation failed for: ", smiles)
      return(tibble::tibble(SMILES = smiles))
    }
    
    json_text <- paste(java_output, collapse = "")
    
    descs <- tryCatch(rjson::fromJSON(json_text),
                      error = function(e) return(NULL))
    
    if (is.null(descs)) {
      message("JSON parse failed for: ", smiles)
      return(tibble::tibble(SMILES = smiles))
    }
    
    descs <- lapply(descs, function(x) if (length(x) == 0) NA else x)
    
    descs_df <- tryCatch(
      as.data.frame(descs, stringsAsFactors = FALSE),
      error = function(e) tibble()
    )
    
    descs_df$SMILES <- smiles
    descs_df
  }
  
  purrr::map_dfr(df$SMILES, calc_by_smiles)
}

               
               

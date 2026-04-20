#' Export aligned feature table from patRoon to NeatMS
#' This function exports the aligned feature table from patRoon to a format compatible with NeatMS. It extracts the necessary information from the XCMSnExp object and the patRoon feature groups, and then combines them into a single dataframe that can be exported as a CSV file for use in NeatMS.
#' @param fGroups The feature groups object from patRoon (default is fGroupsXCMS).
#' @param csv_path The file path where the aligned feature table CSV will be saved (default is "aligned_feature_table.csv").
#' @return A dataframe containing the aligned feature information, which is also saved as a CSV file for NeatMS.
#' @importFrom xcms featureDefinitions featureValues
#' @importFrom dplyr bind_rows left_join select rename filter group_by summarise
#' @importFrom MSnbase sampleNames
#' @importFrom tidyr pivot_longer pivot_wider drop_na
#' @importFrom tibble rownames_to_column
#' @export

export_neatms <- function(fGroups = fGroupsXCMS, csv_path = "aligned_feature_table.csv"){

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
  feature_dataframe <- bind_rows(feature_dataframe, df_filled_peaks)
  
  # Rename the retention time columns of the adjusted rt dataframe
  df_chrom_peaks_fill <- df_chrom_peaks_fill %>%
    rename(rt_adjusted = rt, rtmin_adjusted = rtmin, rtmax_adjusted = rtmax)
  
  # Add the adjusted rt columns to the dataframe containing raw rt
  # To have a dataframe containing both raw and adjusted rt information
  feature_dataframe <- left_join(rownames_to_column(feature_dataframe), rownames_to_column(df_chrom_peaks_fill[,c("rt_adjusted","rtmin_adjusted","rtmax_adjusted")]), by="rowname")
  
  # Remove the rownames as we won't need them
  feature_dataframe$rowname <- NULL
  
  # Retrieve the sample names and store them as a dataframe
  sample_names_df <- as.data.frame(sampleNames(xdata_filled))

#------------------
# Add sample names manually from patRoon data (Drew)
#------------------
  sample_names_df <- as.data.frame(fGroups@features@analysisInfo[["analysis"]])
  sample_names_df$`fGroups@features@analysisInfo[["analysis"]]` <- paste0(sample_names_df$`fGroups@features@analysisInfo[["analysis"]]`, ".mzML")
  
  # Rename the unique column "sample_name"
  colnames(sample_names_df) <- c("sample_name")
  
  # Generate the correct sample ids for matching purposes
  # MSnbase sampleNames() function returns sample names ordered by their ids
  sample_names_df$sample <- seq.int(nrow(sample_names_df))
  
  # Attach the sample names to the main dataframe by matching ids (sample column)
  feature_dataframe <- left_join(feature_dataframe,sample_names_df, by="sample")
  
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
  features_df = rename(features_df, "peak_id"=all_of(peak_colummn_name))

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



#' Import NeatMS results and create YAML file for removal of noisy features
#' This function imports the results from NeatMS, combines them with the metadata from patRoon, and creates a YAML file for the removal of noisy features based on the quality labels assigned by NeatMS.
#' @param ntms_results The file path to the NeatMS results CSV file (default is "ntms_export.csv").
#' @param anaInfo The analysis information dataframe from patRoon.
#' @param feature_dataframe The dataframe containing feature information, including feature IDs and sample names.
#' @param yaml_path The file path where the YAML file for removal of noisy features will be saved (default is "model_session.yml").
#' @return A dataframe summarizing feature quality by group, which is also saved as a CSV file for further analysis.
#' @importFrom dplyr left_join group_by summarise filter distinct select rename
#' @importFrom yaml as.yaml write
#' @importFrom patRoon as.data.table
#' @export
import_neatms <- function(ntms_results = "ntms_export.csv", anaInfo = anaInfo, feature_dataframe = feature_dataframe, yaml_path = "model_session.yml"){

# Import results
neatms_export <- read.csv(ntms_results)

# Select relevant columns from feature_dataframe
feature_dataframe <- feature_dataframe %>%
  select(feature_id, into)

    if (missing(feature_dataframe)) stop("feature_dataframe must be provided")
  if (missing(anaInfo)) stop("anaInfo must be provided")

# Rename and transform columns in NeatMS export
neatms_export <- neatms_export %>%
  rename(
    feature_id = feature.id,
    mz = m.z,
    maxo = height,
    into = area,
    rt = retention.time,
    analysis = sample
  ) %>%
  mutate(rt = rt * 60)  # Convert retention time from minutes to seconds

# Join metadata from anaInfo and feature_dataframe
neatms_export <- neatms_export %>%
  left_join(select(anaInfo, group, analysis), by = "analysis") %>%
  left_join(feature_dataframe, by = "into")

write.csv(neatms_export, "neatms_export_aligned.csv", row.names = FALSE)
  
# Summarize feature quality by group
neatms_export <- neatms_export %>%
  group_by(feature_ID, group) %>%
  summarise(
    quality = any(label == "High_quality"),
    feature = first(feature_id),
    .groups = "drop"
  )



# Filter features with <1 "high-quality"
removeFully <- neatms_export %>%
  group_by(feature) %>%
  filter(sum(quality)/length(quality) == 0) %>%
  distinct(feature) %>%
  select(feature) %>%
  rename(removeFully = feature)

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

featureGroups <- as.data.table(fGroups) %>%
  select(group, ret, mz) %>%
  filter(group %in% removeFully$removeFully)
featureGroups <- list(featureGroups = split(replace(featureGroups, "group", NULL), featureGroups$group))
list <- append(list, featureGroups)

# Convert the list to a YAML object
yaml_object <- as.yaml(list)

# Write the YAML object to a file
write(yaml_object, file = yaml_path)

}
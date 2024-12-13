# ------------------
# GNPS Export
# ------------------

generateGNPSInfo <- function(fGroups, mslists, path) {
  
  # Check inputs for validity
  if (!"xdata" %in% slotNames(fGroups)) stop("fGroups does not contain an xdata slot")
  if (!dir.exists(dirname(path))) stop("Output path does not exist")
  
  # Generate feature table
  featuresDef <- featureDefinitions(fGroups@xdata)
  featuresIntensities <- featureValues(fGroups@xdata, value = "into")
  peakID <- as.character(colnames(fGroups@groups))
  
  # Merge features and intensities, ensuring proper column names
  dataTable <- merge(featuresDef, featuresIntensities, by = 0, all = TRUE)
  colnames(dataTable)[colnames(dataTable) == "Row.names"] <- "Feature_ID"  # Rename merged key column
  dataTable <- dataTable[, !(colnames(dataTable) %in% c("peakidx"))]  # Drop unnecessary columns
  dataTable$peak_ID <- peakID  # Add peak_ID column
  
  # Write feature table to a CSV file
  write.table(dataTable, file.path(path, "output_gnps.csv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Convert mslists to a data.table if needed
  resultsmslists <- as.data.table(mslists)
  
  # Initialize MGF output with metadata
  output_gnps <- vector("list", 2 * length(unique(resultsmslists$group)) + 1)
  output_gnps[[1]] <- paste0("COM=Exported by Drew Szabo on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  index <- 2
  
  # Loop through each unique group
  unique_groups <- unique(resultsmslists$group)
  for (peak_index in seq_along(unique_groups)) {
    peak_ID <- unique_groups[peak_index]
    
    # Extract feature corresponding to the current peak ID
    feature_list <- dataTable[dataTable$peak_ID == peak_ID, ]
    if (nrow(feature_list) == 0) next  # Skip if no matching feature found
    
    # Extract MS1 and MS2 peaks
    MS1peaks <- resultsmslists[group == peak_ID & type == "MS", .(mz, intensity)]
    MS2peaks <- resultsmslists[group == peak_ID & type == "MSMS", .(mz, intensity)]
    
    # Extract metadata for the feature
    ret <- round(feature_list$rtmed, 2)
    mz <- round(feature_list$mzmed, 4)
    numPeaksMS1 <- nrow(MS1peaks)
    numPeaksMS2 <- nrow(MS2peaks)
    feature_ID <- feature_list$Feature_ID
    
    # Generate MS1 entry
    output_gnps[[index]] <- c(
      "BEGIN IONS",
      paste0("FEATURE_ID=", feature_ID),
      paste0("PEAK_ID=", peak_ID),
      "MSLEVEL=1",
      paste0("RTINSECONDS=", ret),
      paste0("PEPMASS=", mz),
      paste0("SCANS=", peak_index),
      paste0("Num peaks=", numPeaksMS1),
      apply(MS1peaks, 1, function(row) paste(round(row["mz"], 6), round(row["intensity"], 0))),
      "END IONS",
      ""
    )
    index <- index + 1
    
    # Generate MS2 entry
    output_gnps[[index]] <- c(
      "BEGIN IONS",
      paste0("FEATURE_ID=", feature_ID),
      paste0("PEAK_ID=", peak_ID),
      "MSLEVEL=2",
      paste0("RTINSECONDS=", ret),
      paste0("PEPMASS=", mz),
      paste0("SCANS=", peak_index),
      paste0("Num peaks=", numPeaksMS2),
      apply(MS2peaks, 1, function(row) paste(round(row["mz"], 6), round(row["intensity"], 0))),
      "END IONS",
      ""
    )
    index <- index + 1
  }
  
  # Write MGF file
  fileConn <- file(file.path(path, "output_gnps.mgf"))
  writeLines(unlist(output_gnps), fileConn)
  close(fileConn)
}

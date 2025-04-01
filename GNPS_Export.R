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
                 paste0("SCANS=", which(dataTable$peak_ID == peak_ID)),
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
                 paste0("SCANS=", which(dataTable$peak_ID == peak_ID)),
                 paste0("Num peaks=", numPeaksMS2),
                 apply(MS2peaks, 1, function(x) paste(round(x[1], 6), round(x[2], 0), sep = " ")),
                 "END IONS", ""), fileConn)
  }
  
  close(fileConn)
}

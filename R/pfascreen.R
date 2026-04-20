#' Export feature list for PFAScreen
#' This function exports a feature list from patRoon in a format compatible with PFAScreen. It extracts the necessary information from the MS lists and feature groups, and then combines them into a single dataframe that can be exported as a CSV file for use in PFAScreen.
#' @param mslists The MS lists object from patRoon (default is mslists).
#' @param path The file path where the output CSV file will be saved (default is "results").
#' @return A dataframe containing the feature information, which is also saved as a CSV file for PFAScreen.
#' @export
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


#' Import PFAScreen results and combine with patRoon metadata
#' This function imports the results from PFAScreen, combines them with the metadata from patRoon, and creates a dataframe that can be used for further analysis. It reads the PFAScreen results from a CSV file, merges them with the feature groups from patRoon, and filters the results based on the number of unique homologues.
#' @param folder The folder path where the PFAScreen results CSV file is located (default is "results/PFAScreen-main"). The function expects a file named "Results_PFAScreen.*\\.csv" within this folder.
#' @return A dataframe containing the combined results from PFAScreen and patRoon, which can be used for further analysis.
#' @export
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
               


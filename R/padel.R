#' Calculate molecular descriptors using PaDEL-Descriptor
#' This function calculates molecular descriptors for a given dataframe containing SMILES strings using the PaDEL-Descriptor Java application. It takes a dataframe with a "SMILES" column, runs the PaDEL-Descriptor tool for each SMILES string, and returns a dataframe with the calculated descriptors. The function also includes error handling to manage cases where the descriptor calculation fails for certain SMILES strings.
#' @param df A dataframe containing a column named "SMILES" with the SMILES strings for which descriptors are to be calculated.
#' @param jar_path The file path to the PaDEL-Descriptor JAR file (default is "descriptor-cli-0.1a-SNAPSHOT-all.jar").
#' @return A dataframe containing the calculated molecular descriptors for each SMILES string, along with the original SMILES string.
#' @export
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
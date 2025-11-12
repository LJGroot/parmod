add_dummies_help <- function(){
  cat(r"(
## FUNCTION `add_dummies()`
# This function can be used to add dummy variables to you IPD data set. 
  
# It assumes that your IPD is stored as raw data in long format with a study identifier (preferably labeled as "study").
# The function has the following arguments:
# df = provide data frame
# study_col = provide column that holds study ID, default is "study"
# ref = reference group, first study by default
# prefix = preferred prefix for names of dummy variables, defaults to "d"
# name_by_value = use numbers 1 through k as dummy suffix (TRUE), 
#                 or use original study ID as suffix (FALSE)
  )")
  }

add_dummies <- function(df, 
                        study_col = "study", 
                        ref = NULL, 
                        prefix = "d", 
                        name_by_value = FALSE) {
  # error for wrong column name
  if (!study_col %in% names(df)) {
    stop("Column '",
study_col, 
"' not found in data frame. 
Run add_dummies_help() for help on using this function")
  }
  # full-length vector of values (not unique!)
  studies_vec <- as.character(df[[study_col]])
  # unique values as vector
  unique_studies <- unique(studies_vec)
    # error for when column has no values
  if (length(unique_studies) == 0) {
    stop("Study column contains no values. Run add_dummies_help() for help on using this function")
  }
  # if ref argument is null, take first unique entry as ref
  # if ref is provided but not present, give error
  if (is.null(ref)) {
    ref <- unique_studies[1]
  } else if (!ref %in% unique_studies) {
    stop("Reference '", ref, "' not found among study values:",
         paste(head(unique_studies, 10), collapse = ", "))
  }
  # take all levels minus reference group
  levels_to_dummy <- unique_studies[unique_studies != ref]
  # throw error if remaining levels is zero (only 1 unique study in data)
  if (length(levels_to_dummy) == 0) {
    message("only the reference level '", ref, "' present - no dummy columns created.")
    return(df)
  }
  # generate matrix of dummies
  dummy_mat <- vapply(levels_to_dummy,
                      FUN = function(val) as.integer(studies_vec == val),
                      FUN.VALUE = integer(length(studies_vec)),
                      USE.NAMES = FALSE)
  # convert to data.frame
  dummy_df <- as.data.frame(dummy_mat, stringsAsFactors = FALSE)
  # give names to dummies using provided of default prefix, and using
  # numbers or original study identifiers
  if (name_by_value) {
    colnames(dummy_df) <- paste0(prefix, make.names(levels_to_dummy))
  } else {
    colnames(dummy_df) <- paste0(prefix, seq_along(levels_to_dummy))
  }
  
  out <- df
  out[names(dummy_df)] <- dummy_df
  out
}

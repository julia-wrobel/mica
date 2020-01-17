#' Vectorize and name nifti object
#'
#' Function used to read in nifti object, vectorize, and place in a
#' dataframe. Option to whitestripe normalize image if desired.
#'
#' @param filepath Path to where nifti object is stored.
#' @param site_scan_id Character valued ID for site, scan number, and id in the format
#' site_scannumber_id.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom WhiteStripe whitestripe whitestripe_norm
#' @importFrom neurobase readnii
#' @importFrom dplyr filter
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

vectorize_image <- function(filepath, site_scan_id = NULL,
                            filter_skull = TRUE, ...){
  nifti_object = readnii(filepath)

  df = data.frame(value = c(nifti_object))
  colnames(df) = "intensity"

  df$id = strsplit(site_scan_id, "_")[[1]][3]
  df$site = strsplit(site_scan_id, "_")[[1]][1]
  df$scan = strsplit(site_scan_id, "_")[[1]][2]
  df$voxel_position = row.names(df)

  filter(df, intensity > 0)
  #filter(df, intensity != round(0, 1)) # 0 is the value for skull
}

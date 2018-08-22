#' Vectorize and name nifti object
#'
#' Function used to read in nifti object, vectorize, and place in a
#' dataframe. Option to whitestripe normalize image if desired.
#'
#' @param filepath Path to where nifti object is stored.
#' @param site_scan_id Character valued ID for site, scan number, and id in the format
#' site_scannumber_id.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized.
#' @param type If \code{white_stripe = TRUE}, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
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

vectorize_image <- function(filepath, site_scan_id = NULL, white_stripe = FALSE, type = NULL, ...){
  nifti_object = readnii(filepath)

  if(white_stripe){
    ind = whitestripe(img = nifti_object, type = type,
                      stripped = TRUE)$whitestripe.ind
    nifti_object = whitestripe_norm(nifti_object, indices = ind)
  }

  df = data.frame(value = c(nifti_object))
  colnames(df) = "intensity"

  df$id = strsplit(site_scan_id, "_")[[1]][3]
  df$site = strsplit(site_scan_id, "_")[[1]][1]
  df$scan = strsplit(site_scan_id, "_")[[1]][2]

  filter(df, intensity > min(intensity))
}

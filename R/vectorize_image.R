#' vectorize and name nifti object
#'
#' Function used to read in nifti object, vectorize, and place in a
#' dataframe. Option to whitestripe normalize image if desired.
#'
#' @param filepath Path to where nifti object is stored.
#' @param site_scan Character valued ID for site and scan number in the format site_scannumber.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized.
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

vectorize_image <- function(filepath, site_scan = NULL, white_stripe = FALSE, ...){
  nifti_object = readnii(filepath)

  if(white_stripe){
    ind = whitestripe(img = nifti_object, type = "T1",
                      stripped = TRUE)$whitestripe.ind
    nifti_object = whitestripe_norm(nifti_object, indices = ind)
  }

  df = data.frame(value = c(nifti_object))

  # # probably need to edit this part
  # if(!is.null(name_indices)){
  #   sitename = substring(filepath, name_indices[1], name_indices[1]+2)
  #   scan_no = substring(str_split(filepath, "Scan_")[[1]][2], 1, 1)
  #
  #   colnames(df) = paste0(sitename, "_", scan_no)
  # }

  colnames(df) = "intensity"

  df$id = strsplit(site_scan, "_")[[1]][3]
  df$site = strsplit(site_scan, "_")[[1]][1]
  df$scan = strsplit(site_scan, "_")[[1]][2]

  filter(df, intensity > 0)
}

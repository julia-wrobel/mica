#' Normalize and save nifti object
#'
#' Function used to mica-normalize nifti object based on inverse warping functions.
#'
#' @param filepath Path to where nifti object is stored.
#' @param subject Unique identifiers for the nifti object.
#' @param scan_id Unique identifiers for the nifti object.
#' @param data Data frame that includes the variable \code{h_inv}, which is the
#' inverse warping function for normalization.
#' @param outpath Directory where mica normalized nifti object will be stored.
#' @param white_striped Has the data been intensity normalized using White Stripe?
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom neurobase readnii niftiarr writenii
#'
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

normalize_image = function(filepath, subject, scan_id, data, outpath,
                           white_striped, ...){
  nifti_unnorm = readnii(filepath)

  nifti_vec = data.frame(intensity = c(nifti_unnorm))
  # getmode <- function(v) {
  #   uniqv <- unique(v)
  #   uniqv[which.max(tabulate(match(v, uniqv)))]
  # }

  # don't want to recalculate mode because it is computationally intensive
  mode = unique(data$mode)
  nifti_vec$intensity[which(nifti_vec$intensity != mode)] = data[["h_inv"]]

  nifti_norm = niftiarr(img = nifti_unnorm, nifti_vec$intensity)
  new_filename = paste0(outpath, "/", subject, "_", scan_id, "_mica.nii.gz")
  writenii(nifti_norm, filename = new_filename)
}

#' Normalize and save nifti object
#'
#' Function used to mica-normalize nifti object based on inverse warping functions.
#'
#' @param filepath Path to where nifti object is stored.
#' @param scan_id Unique identifiers for the nifti objecgt.
#' @param data Data frame that includes the variable \code{h_inv}, which is the
#' inverse warping function for normalization.
#' @param outpath Directory where mica normalized nifti object will be stored.
#' @param white_striped Has the data been intensity normalized using White Stripe?
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom neurobase readnii niftiarr writenii
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

# need to convert back whitestripe images... are those = 0 in the same places? ignore for now
normalize_image = function(filepath, scan_id, data, outpath,
                           white_striped, ...){
  nifti_unnorm = readnii(filepath)

  nifti_vec = data.frame(intensity = c(nifti_unnorm))

  if(white_striped){
    nifti_vec$intensity[which(nifti_vec$intensity != min(nifti_vec$intensity))] = data[["h_inv"]]
  }else{
    nifti_vec$intensity[which(nifti_vec$intensity != 0)] = data[["h_inv"]]
  }

  nifti_norm = niftiarr(img = nifti_unnorm, nifti_vec$intensity)
  new_filename = paste0(outpath, "/", scan_id, "_mica.nii.gz")
  writenii(nifti_norm, filename = new_filename)
}

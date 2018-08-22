#' Normalize and save nifti object
#'
#' Function used to mica-normalize nifti object based on inverse warping functions.
#'
#' @param filepath Path to where nifti object is stored.
#' @param id Unique identifiers for the nifti object.
#' @param data Data frame that includes the variable \code{h_inv}, which is the
#' inverse warping function for normalization.
#' @param outpath Directory where mica normalized nifti object will be stored.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom neurobase readnii niftiarr writenii
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

# need to convert back whitestripe images... are those = 0 in the same places? ignore for now
normalize_image = function(filepath, id, data, outpath, ...){

  #return(list(filepath = filepath, id = id, data = data))

  nifti_unnorm = readnii(filepath)

  nifti_vec = data.frame(intensity = c(nifti_unnorm))
  nifti_vec$intensity[which(nifti_vec$intensity > 0)] = data[["h_inv"]]

  nifti_norm = niftiarr(img = nifti_unnorm, nifti_vec$intensity)
  new_filename = paste0(outpath, "/", id, "_mica.nii.gz")
  writenii(nifti_norm, filename = new_filename)
}

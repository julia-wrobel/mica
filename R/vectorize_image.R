#' Vectorize and name nifti object
#'
#' Function used to read in nifti object, vectorize, and place in a
#' dataframe. Option to whitestripe normalize image if desired.
#'
#' @param filepath Path to where nifti object is stored.
#' @param subj_scan_scanner Character valued ID for subject, scan number, and scanner, in the format
#' subj_scan_scanner.
#' @param white_striped Has the data been intensity normalized using White Stripe? Defaults to FALSE.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom neurobase readnii
#' @importFrom dplyr filter
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

vectorize_image <- function(filepath, subj_scan_scanner,
                            white_striped, ...){
  nifti_object = readnii(filepath)

  df = data.frame(intensity = c(nifti_object))

  df$subject = strsplit(subj_scan_scanner, "_")[[1]][1]
  df$scan_id = strsplit(subj_scan_scanner, "_")[[1]][2]
  df$scanner = strsplit(subj_scan_scanner, "_")[[1]][3]
  df$voxel_position = row.names(df)

  # want it to filter the zeros, add minimum value, store minimum value
  if(white_striped){
    df = filter(df, intensity != min(intensity))
  }else if(min(df$intensity) < 0){
    df = filter(df, intensity != 0)
    df = mutate(df, intensity = intensity - min(intensity))
  }else{
    df = filter(df, intensity != 0)
  }

  df
}


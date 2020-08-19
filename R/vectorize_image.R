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
#' @author Julia Wrobel \email{julia.wrobel@@cuanschutz.edu}
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

  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  # want it to filter the zeros, add minimum value, store minimum value
  mode = getmode(df$intensity)
  df = mutate(df, mode = mode)
  df = filter(df, intensity != mode)

  df
}


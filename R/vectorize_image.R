#' Vectorize and name nifti object
#'
#' Function used to read in nifti object, vectorize, and place in a
#' dataframe. Option to whitestripe normalize image if desired.
#'
#' @param filepath Path to where nifti object is stored.
#' @param site_id Character valued ID for site, and id in the format
#' site_id.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom neurobase readnii
#' @importFrom dplyr filter
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @return a data frame with a single vectorized image.
#' @export

vectorize_image <- function(filepath, site_id,
                            filter_skull = TRUE, ...){
  nifti_object = readnii(filepath)

  df = data.frame(intensity = c(nifti_object))

  df$id = strsplit(site_id, "_")[[1]][2]
  df$site = strsplit(site_id, "_")[[1]][1]
  df$voxel_position = row.names(df)

  df = filter(df, intensity != 0)
  # want it to filter the zeros, add minimym value, store minimum value

  if(min(df$intensity) < 0){
    df = mutate(df, intensity = intensity - min(intensity))
  }

  df
}


#' create labeled dataframe of vectorized voxel intensities for a list of nifti objects
#'
#' Function used along with \code{mica::vectorize_image()} to vectorize images and place in a
#' dataframe, where names for each image are drawn from the filepath of the
#' image.
#'
#' @param filepaths List of paths to where nifti objects are stored.
#' @param sitenames List of character valued IDs for imaging site.
#' @param scan_nos List of character valued IDs for scan number. If NULL, scan number will be assumed to be 1.
#' @param ids List of unique identifiers for each image.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom purrr map2 map_df
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @examples
##' \dontrun{
##' filenames = list.files(path = "../Documents/", pattern = "flair_reg", recursive = TRUE)
##' filenames = paste0("../Documents/", filenames)
##' sitenames = lapply(filenames, function(f) substring(f, 64, 66))
##' scan_nos = lapply(filenames, function(f) substring(str_split(f, "Scan_")[[1]][2], 1, 1))
##' ids = 1:14
##'
##' intensity_data = make_intensities_df(filenames, sitenames, scan_nos, ids)
##' }
#'
#' @return a data frame of vectorized images with columns \code{intensity}, \code{id},
#' \code{site}, and \code{scan}, \code{cdf}.
#' @export

make_intensity_df <- function(filepaths, sitenames, scan_nos, ids, ...){

  # for testing filepaths, sitenames, scan_nos, ids should all have the same length and should be lists
  #
  site_scan_id = paste(sitenames, scan_nos, ids, sep = "_")
  intensities = map2(filepaths, site_scan_id, vectorize_image)

  intensities = map_df(intensities, rbind)
  #mutate(intensities, cdf = )
}

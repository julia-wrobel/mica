#' create labeled dataframe of vectorized voxel intensities for a list of nifti objects
#'
#' Function used along with \code{mica::vectorize_image()} to vectorize images and place in a
#' dataframe, where names for each image are drawn from the filepath of the
#' image.
#'
#' @param filepaths List of paths to where nifti objects are stored.
#' @param ids List of unique identifiers for each image.
#' @param sitenames List of character valued IDs for imaging site.
#' @param scan_nos List of character valued IDs for scan number. If NULL, scan number will
#' be assumed to be 1.
#' @param white_stripe If \code{TRUE} image will be white stripe normalized
#' before it is vectorized
#' @param type If \code{white_stripe = TRUE}, user must specify the type of image, from options
#' \code{type = c("T1", "T2", "FA", "MD", "first", "last", "largest")}.
#' @param ... Additional arguments passed to or from other functions.
#'
#' @importFrom purrr map2 map_df
#' @importFrom tibble as_tibble
#'
#' @author Julia Wrobel \email{jw3134@@cumc.columbia.edu}
#'
#' @examples
##' \dontrun{
##' filenames = list.files(path = "../Documents/", pattern = "flair_reg", recursive = TRUE)
##' filenames = paste0("../Documents/", filenames)
##' sitenames = lapply(filenames, function(f) substring(f, 64, 66))
##' scan_nos = lapply(filenames, function(f) substring(str_split(f, "Scan_")[[1]][2], 1, 1))
##' ids = paste(sitenames, scan_nos, sep = "_")
##'
##' intensity_data = make_intensities_df(filenames, ids, sitenames, scan_nos)
##' }
#'
#' @return a data frame of vectorized images with columns \code{intensity}, \code{id},
#' \code{site}, and \code{scan}, \code{cdf}.
#' @export

make_intensity_df <- function(filepaths, ids, sitenames = NULL, scan_nos = NULL,
                              white_stripe = FALSE, type = NULL, ...){

  # for testing filepaths, sitenames, scan_nos, ids should all have the same length and
    #should be lists
  # also test that if sitenames and scan_nos is null you don't through an error
  if(is.null(sitenames)) sitenames = rep("site", length(filepaths))
  if(is.null(scan_nos)) scan_nos = rep("1", length(filepaths))

  site_scan_id = paste(sitenames, scan_nos, ids, sep = "_")
  intensities = map2(filepaths, site_scan_id, vectorize_image,
                     white_stripe = white_stripe, type = type)

  intensities = map_df(intensities, rbind)

  as_tibble(intensities)
}

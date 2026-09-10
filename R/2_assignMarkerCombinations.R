#' Assing each cell to a combination of gated markers.
#'
#' @param data_markerPos A data.frame where each row is a cell and each column is a gated marker describing
#' whether a cell is included in the gate, as constructed by the \code{\link{getMarkerPositivity}} function.
#'
#' @param mode Character string specifying the level of detail in the polyfunctionality assignment.
#' If \code{simple}, the output includes only \code{MarkerComb} (exact marker combinations),
#' \code{MarkerCount} (number of functions), and \code{Functionality} (mono/polyfunctional).
#' If \code{exhaustive}, the output additionally contains all "atleast_n" columns, which are logical
#' describing whether cells have at least 1, 2, 3, ... \eqn{n} markers in any combination. Default = \code{simple}.
#'
#' @returns The input data.frame with three extra columns describing
#'     \itemize{
#'       \item{\code{MarkerComb}: A factor describing which marker combination a specific cell has. Describes
#'       the exact combinations, which are mutually exclusive (e.g., 4 markers → 2^4 = 16 combinations).}
#'       \item{\code{MarkerCount}: An integer describing the number of functions of the cell (count of positive markers).}
#'       \item{\code{Functionality}: A factor grouping \code{MarkerCount)} into three categories; No cytokines (\code{MarkerCount == 0}),
#'       Monofunctional (\code{MarkerCount == 1}), or Polyfunctional (\code{MarkerCount > 1}).}
#'     }
#'
#' @seealso \code{\link{getMarkerPositivity}}, \code{\link{calcPolyfunctionality}}
#' @importFrom utils combn
#' @importFrom rlang .data :=
#' @examples
#' # Define path to data
#' path_data <- system.file("extdata", package = "polyICSFlow")
#'
#' # Read fcs files as flowSet
#' fs <- flowCore::read.flowSet(path = path_data,
#'                              pattern = ".fcs",
#'                              truncate_max_range = FALSE,
#'                              transformation = FALSE)
#'
#' # Create GatingSet
#' gs <- flowWorkspace::GatingSet(fs)
#'
#' # Apply GatingTemplate
#' gt <- openCyto::gatingTemplate(file.path(path_data, "gatingTemp_cytokines.csv"))
#' openCyto::gt_gating(x = gt, y = gs)
#'
#' # Get marker positivity per cell
#' df_markerPos <- getMarkerPositivity(x = gs,
#'                                     gate_names = c("IFN+","TNFa+","IL2+","CD107a+"))
#'
#' # Assign marker combinations to each cell based on marker positivity
#' # Simple mode
#' df_markerComb <- assignMarkerCombinations(data_markerPos = df_markerPos,
#'                                           mode = "simple")
#' # Exhaustive mode
#' df_markerComb <- assignMarkerCombinations(data_markerPos = df_markerPos,
#'                                           mode = "exhaustive")
#' @export
assignMarkerCombinations <- function(data_markerPos, mode = c("simple","exhaustive")){

  mode <- match.arg(arg = mode, choices = c("simple", "exhaustive"))

  markers <- names(data_markerPos)
  n_markers <- length(markers)
  # get order of combinations
  unique_combinations <- .findUniqueMarkerCombinations(marker_names = markers, full_names = TRUE)

  msg <- paste0("Assigning cells to marker combinations (", length(unique_combinations), " total from ", n_markers, " markers)...")
  message(msg)

  df_simple <- data_markerPos %>%
    dplyr::mutate(MarkerComb = apply(dplyr::select(data_markerPos, dplyr::where(is.logical)),
                                     1,
                                     function(x) {paste0(names(x),
                                                         ifelse(x, "+", "-"),
                                                         collapse = "")}),
                  MarkerCount = as.integer(rowSums(dplyr::across(dplyr::where(is.logical)))),
                  Functionality = dplyr::case_when(MarkerCount == 0 ~ "None",
                                                   MarkerCount == 1 ~ "Monofunctional",
                                                   MarkerCount > 1 ~ "Polyfunctional"))%>%
    dplyr::mutate(MarkerComb = factor(.data$MarkerComb,levels = unique_combinations),
                  Functionality = factor(.data$Functionality, levels = c("None","Monofunctional","Polyfunctional")))

  message("Done!")

  if(mode == "simple"){

    df_final <- df_simple

  }else if(mode == "exhaustive"){

    df_exhaustive <- df_simple %>%
      .getAtLeast_n() %>%
      .getPolyfunctionalAtLeast_marker()

    df_final <- df_exhaustive

    message("Done!")

  }


  return(df_final)

}

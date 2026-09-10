#' Plot a heatmap of all possible marker combinations.
#'
#' @param marker_names A character vector specifying the marker names used to
#'   generate all possible combinations.
#' @param palette_count A named character vector mapping the number of markers
#'   per cell to colors. If \code{NULL} (default), a qualitative palette is generated with
#'   \code{\link[scales]{brewer_pal}}.
#' @param palette_comb A named character vector mapping marker combination
#'   names to colors. If \code{NULL} (default), a qualitative palette is generated with
#'   \code{CATALYST:::.cluster_cols}.
#' @param orientation A string specifying the orientation of the heatmap, either vertical or horizontal.
#'  Default = \code{"vertical"}.
#' @param full_names Logical. If TRUE (default) the marker combinations are displayed with their full
#'  names (i.e., negative markers are also present). If FALSE, only the positive markers are displayed.
#'
#' @return A \code{\link[ComplexHeatmap]{Heatmap}} object visualizing all marker
#'   combinations. The heatmap includes two layers of annotation:
#'   \itemize{
#'     \item \code{MarkerCount}: A marker count annotation (colored
#'       with \code{palette_count}) showing the number of markers in each combination.
#'     \item \code{MarkerComb}: A combination annotation (colored
#'       with \code{palette_comb}) showing the exact subset of markers present.
#'   }
#'
#' @seealso \code{\link[ComplexHeatmap]{Heatmap}}
#' @export
#' @examples
#' # Vertical heatmap
#' plotMarkerCombHeatmap(marker_names = c("IFN","TNFa","IL2","CD107a"),
#'                       orientation = "vertical")
#'
#' # Horizontal heatmap
#' plotMarkerCombHeatmap(marker_names = c("IFN","TNFa","IL2","CD107a"),
#'                       orientation = "horizontal")
#' @export
plotMarkerCombHeatmap <- function(marker_names, orientation = c("vertical","horizontal"), palette_count = NULL, palette_comb = NULL, full_names = TRUE){

  orientation <- match.arg(orientation)

  levels_combinations <- .findUniqueMarkerCombinations(marker_names, full_names = full_names)

  # get matrix of TRUE/FALSE
  matrix_legend <- outer(marker_names,
                         levels_combinations,
                         Vectorize(function(marker, combo) {
                           grepl(paste0(marker, "\\+"), combo)
                         }))
  # Convert to 0/1
  matrix_legend <- 1L * matrix_legend

  # Create column annotation
  df_annotation <- data.frame(MarkerComb = factor(levels_combinations),
                              MarkerCount = factor(colSums(matrix_legend)))
  unique_grouped <- unique(df_annotation$MarkerCount)
  unique_individual <- unique(df_annotation$MarkerComb)

  # Add col and row names
  colnames(matrix_legend) <- levels_combinations
  rownames(matrix_legend) <- marker_names

  # Define palette
  if(is.null(palette_count)){
    palette_count <- c("grey85",scales::brewer_pal(type = "qual")(length(unique_grouped)-1))
  }
  if(is.null(palette_comb)){
    palette_comb <- .getNamedPaletteMarkerComb(marker_names = marker_names, full_names = full_names)
  }
  names(palette_count) <- unique_grouped

  if(orientation == "horizontal"){

    annot <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation,
                                               col = list(MarkerCount = palette_count,
                                                          MarkerComb = palette_comb),
                                               show_legend = c(FALSE,FALSE),
                                               show_annotation_name = c(TRUE,TRUE),
                                               annotation_name_gp = grid::gpar(fontsize = 8),
                                               which = "column")
    p <- ComplexHeatmap::Heatmap(matrix_legend,
                                 show_row_names = TRUE,
                                 row_names_side = "left",
                                 top_annotation = annot,
                                 cluster_rows = FALSE,
                                 col = c("0" = "grey96", "1" = "black"),
                                 cluster_columns = FALSE,
                                 show_heatmap_legend = FALSE,
                                 column_split = df_annotation$MarkerCount,
                                 column_title = NULL,
                                 row_title = NULL,
                                 show_column_names = TRUE,
                                 rect_gp = grid::gpar(col = "white", lwd = 1),
                                 row_names_gp = grid::gpar(fontsize = 8),
                                 column_names_gp = grid::gpar(fontsize = 8))


  }else if(orientation == "vertical"){

    annot <- ComplexHeatmap::HeatmapAnnotation(df = df_annotation,
                                               col = list(MarkerCount = palette_count,
                                                          MarkerComb = palette_comb),
                                               show_legend = FALSE,
                                               show_annotation_name = TRUE,
                                               annotation_name_gp = grid::gpar(fontsize = 8),
                                               annotation_name_align = TRUE,
                                               which = "row")

    p <- ComplexHeatmap::Heatmap(t(matrix_legend),
                                 show_row_names = TRUE,
                                 row_names_side = "left",
                                 right_annotation = annot,
                                 cluster_rows = FALSE,
                                 col = c("0" = "grey96", "1" = "black"),
                                 cluster_columns = FALSE,
                                 show_heatmap_legend = FALSE,
                                 row_split = df_annotation$MarkerCount,
                                 column_title = NULL,
                                 row_title = NULL,
                                 show_column_names = TRUE,
                                 column_names_side = "top",
                                 column_names_rot = 45,
                                 column_names_centered = FALSE,
                                 rect_gp = grid::gpar(col = "white", lwd = 1),
                                 row_names_gp = grid::gpar(fontsize = 8),
                                 column_names_gp = grid::gpar(fontsize = 8))
  }

  return(p)

}
#'
#' Plot counts of cells expressing exactly 1, 2, ..., \eqn{n} markers.
#'
#' Plot a stacked bar plot of cell counts within each marker combination for cells positive for exactly 1, 2, ..., \eqn{n} markers.
#'
#' @param data_markerComb A dataframe generated by the \code{\link{assignMarkerCombinations}} function,
#'   containing marker combination data.
#' @param cell_subset_params List with parameters describing what subset of cells (cellTypes_plot) should be plotted,
#'  from a vector of cell types (cellTypes). Default = \code{list(cellTypes = NULL, cellTypes_plot = NULL)},
#'  i.e. all cells in the dataset are plotted.
#' @param palette_comb A named character vector mapping marker combination
#'   names (as produced by \code{\link{assignMarkerCombinations}}) to colors.
#'   If \code{NULL} (default), a qualitative palette is generated with \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'   and \code{CATALYST:::.cluster_cols}.
#'
#' @returns A ggplot.
#' @export
#'
#' @seealso \code{\link{plotCellCountsAtLeast}}
#'
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
#' df_markerComb <- assignMarkerCombinations(data_markerPos = df_markerPos,
#'                                           mode = "simple")
#'
#' # Plot all cells
#' plotCellCountsExactly(data_markerComb = df_markerComb)
#'
#' # Plot only subset of data (metacluster 2 and 4)
#' metacluster_labels <- readRDS(file.path(path_data, "cell_labels.rds"))
#' plotCellCountsExactly(data_markerComb = df_markerComb,
#'                       cell_subset_params = list(cellTypes = metacluster_labels,
#'                                                 cellTypes_plot = c(2,4)))
plotCellCountsExactly <- function(data_markerComb, cell_subset_params = list(cellTypes = NULL,cellTypes_plot = NULL), palette_comb = NULL){

  markers <- .getMarkerNames(data_markerComb)

  subset_result <- .subsetDataCellTypes(data_markerComb,cell_subset_params = cell_subset_params)
  df <- subset_result$dataframe
  plot_caption <- subset_result$plot_caption

  df_plot <- df %>%
    dplyr::filter(.data$MarkerCount != 0)%>%
    dplyr::select(-c(!!!rlang::syms(markers),.data$Functionality))%>%
    dplyr::group_by(.data$MarkerComb,.data$MarkerCount)%>%
    dplyr::count()%>%
    dplyr::ungroup() %>%
    tidyr::pivot_longer(cols = c(.data$MarkerCount),
                        names_to = "facet")

  if(is.null(palette_comb)){

    palette_comb <- .getNamedPaletteMarkerComb(markers)

  }
  Combs_present <- intersect(names(palette_comb),unique(df_plot$MarkerComb))
  palette_comb <- palette_comb[Combs_present]

  p <- df_plot %>%
    ggplot(aes(x = .data$value,
               y = .data$n,
               fill = .data$MarkerComb))+
    ggplot2::geom_bar(stat = "identity")+
    ggplot2::facet_grid(~facet,
                        labeller = ggplot2::as_labeller(c("MarkerCount" = "Exactly \n n markers")),
                        scales = "free_x",
                        space = "free_x")+
    ggplot2::scale_fill_manual(values = palette_comb,
                      labels = .markerCombtoMarkdown(names(palette_comb)),
                      name = "Marker combination")+
    ggplot2::theme_bw()+
    theme(legend.text = ggtext::element_markdown(hjust = .5,
                                                 face = "bold"),
          strip.text.x = element_text(face = "bold"),
          plot.caption = element_text(size = 8))+
    ylab("Cell count")+
    xlab(NULL)+
    labs(caption = plot_caption)

  return(p)

}


#' Plot counts of cells expressing at least 1, 2, ..., \eqn{n} markers.
#'
#' Plot a stacked bar plot of cell counts by marker combination, faceted by the minimum number of markers
#' co-expressed (\eqn{\geq 1, \geq 2, \ldots, \geq n}). Within each facet, the x-axis shows all marker or marker-combination categories
#' corresponding to that level (e.g., single markers for \eqn{>=1}, two-marker combinations for \eqn{>=2}, etc.).
#'
#' @param data_markerComb A dataframe generated by the \code{\link{assignMarkerCombinations}} function,
#'   containing marker combination data.
#' @param cell_subset_params List with parameters describing what subset of cells (cellTypes_plot) should be plotted,
#'  from a vector of cell types (cellTypes). Default = \code{list(cellTypes = NULL, cellTypes_plot = NULL)},
#'  i.e. all cells in the dataset are plotted.
#' @param palette_comb A named character vector mapping marker combination
#'   names (as produced by \code{\link{assignMarkerCombinations}}) to colors.
#'   If \code{NULL} (default), a qualitative palette is generated with \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'   and \code{CATALYST:::.cluster_cols}.
#'
#' @returns A ggplot.
#' @export
#' @seealso \code{\link{plotCellCountsExactly}}
#' @importFrom stats hclust
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
#' df_markerComb <- assignMarkerCombinations(data_markerPos = df_markerPos,
#'                                           mode = "simple")
#'
#' # Plot all cells
#' plotCellCountsAtLeast(data_markerComb = df_markerComb)
#'
#' # Plot only subset of data (metacluster 2 and 4)
#' metacluster_labels <- readRDS(file.path(path_data, "cell_labels.rds"))
#' plotCellCountsAtLeast(data_markerComb = df_markerComb,
#'                       cell_subset_params = list(cellTypes = metacluster_labels,
#'                                                 cellTypes_plot = c(2,4)))
plotCellCountsAtLeast <- function(data_markerComb, cell_subset_params = list(cellTypes = NULL,cellTypes_plot = NULL), palette_comb = NULL){

  markers <- .getMarkerNames(data_markerComb)
  n_markers <- length(markers)

  # if mode = "simple" was used, i.e. at least columns missing, add them
  if(!any(grepl("atleast",names(data_markerComb)))){
    data_markerComb <- .getAtLeast_n(data_markerComb)
    message("Done!")
  }

  subset_result <- .subsetDataCellTypes(data_markerComb,cell_subset_params = cell_subset_params)
  df <- subset_result$dataframe
  plot_caption <- subset_result$plot_caption

  df_plot <- df %>%
    dplyr::select(c(.data$MarkerComb,dplyr::starts_with("atleast_")))%>%
    tidyr::pivot_longer(cols = -.data$MarkerComb,
                        names_to = "Marker")%>%
    dplyr::filter(.data$value)%>%
    dplyr::group_by(.data$Marker,.data$MarkerComb)%>%
    dplyr::count()%>%
    dplyr::ungroup() %>%
    dplyr::mutate(AtLeast = stringr::str_split_i(.data$Marker, "_", 2),
                  xlab = stringr::str_split_i(.data$Marker, "_", 3))%>%
    dplyr::mutate(xlab = forcats::fct_reorder(.data$xlab, .data$n, .fun = sum, .desc = TRUE))

  if(is.null(palette_comb)){

    palette_comb <- .getNamedPaletteMarkerComb(marker_names = markers)
  }
  Combs_present <- intersect(names(palette_comb),unique(df_plot$MarkerComb))
  palette_comb <- palette_comb[Combs_present]

  p <- df_plot %>%
    ggplot(aes(x = .data$xlab,
               y = .data$n,
               fill = .data$MarkerComb))+
    ggplot2::geom_bar(stat = "identity")+
    ggplot2::facet_grid(~AtLeast,
                        labeller = ggplot2::labeller(AtLeast = function(x) paste("At least \n", x, ifelse(x > 1, "markers","marker"))),
                        scales = "free_x",
                        space = "free_x")+
    ggplot2::scale_fill_manual(values = palette_comb,
                      labels = .markerCombtoMarkdown(names(palette_comb)),
                      name = "Marker combination")+
    ggplot2::theme_bw()+
    theme(legend.text = ggtext::element_markdown(hjust = .5,
                                                 face = "bold"),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = .5),
          strip.text.x = element_text(face = "bold"),
          plot.caption = element_text(size = 8))+
    ylab("Cell count")+
    xlab(NULL)+
    labs(caption = plot_caption)

  return(p)

}


#' Plot 2D scatterplots of all marker combinations.
#'
#' Plots 2D scatter plots of all possible marker combinations, where each panel (facet) corresponds to one of the
#' \eqn{2^n} possible combinations of \eqn{n} markers.
#' The top-left panel represents cells negative for all markers, while the bottom-right represents cells positive for all markers,
#' with intermediate panels showing mixed combinations.
#' Within each panel, the x-axis denotes the \eqn{n} markers and the y-axis shows the transformed MFI values.
#' Each point corresponds to a single cell, with cells belonging to the combination represented in that facet
#' highlighted in color, and all other cells displayed in grey background for reference.
#'
#' @param input Either a matrix of marker expression or an object of class
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}}, \code{\link[flowWorkspace]{GatingSet}},
#'  or \code{\link[flowCore]{flowSet}}, from which the expression data can be extracted.
#'  If input is either GatingSet or flowSet, an aggregate of the data is created with the
#'  \code{\link[FlowSOM]{AggregateFlowFrames}} function.
#' @param assigned_MarkerCombs A vector of assigned marker combinations per
#'   cell in \code{x}, as created by \code{\link{assignMarkerCombinations}}.
#' @param cell_subset_params List with parameters describing what subset of cells (cellTypes_plot) should be plotted,
#'  from a vector of cell types (cellTypes). Default = \code{list(cellTypes = NULL, cellTypes_plot = NULL)},
#'  i.e. all cells in the dataset are plotted.
#' @param palette_comb A named character vector mapping marker combination
#'   names (as produced by \code{\link{assignMarkerCombinations}}) to colors.
#'   If \code{NULL} (default), a qualitative palette is generated with \code{\link[grDevices:colorRampPalette]{colorRampPalette}}
#'   and \code{CATALYST:::.cluster_cols}.
#' @param return_list Logical. If FALSE (default), the plots are combined by cowplot. If TRUE, a list of
#'   plots are returned to allow further customization.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object visualizing 2D scatterplots of all
#'   marker combinations.
#'
#' @importFrom ggplot2 ggplot aes theme ylab xlab element_line element_text
#'   scale_color_manual scale_x_continuous geom_point expansion ggtitle labs theme_classic
#' @importFrom stats dnorm rnorm runif
#' @importFrom methods is
#' @importFrom SummarizedExperiment assay
#' @export
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
#' df_markerComb <- assignMarkerCombinations(data_markerPos = df_markerPos)
#'
#' # Plot 2D scatter of marker combinations in dataset (all data)
#' # From gatingSet
#' plotMarkerComb2DScatters(input = gs,
#'                          assigned_MarkerCombs = df_markerComb$MarkerComb)
#'
#' # From flowSet
#' plotMarkerComb2DScatters(input = fs,
#'                          assigned_MarkerCombs = df_markerComb$MarkerComb)
#'
#' # Plot only subset of data (metacluster 2 and 4)
#' metacluster_labels <- readRDS(file.path(path_data, "cell_labels.rds"))
#' plotMarkerComb2DScatters(input = fs,
#'                          assigned_MarkerCombs = df_markerComb$MarkerComb,
#'                          cell_subset_params = list(cellTypes = metacluster_labels,
#'                                                    cellTypes_plot = c(2,4)))
#' # As list
#' plot_list <- plotMarkerComb2DScatters(input = fs,
#'                                       assigned_MarkerCombs = df_markerComb$MarkerComb,
#'                                       return_list = TRUE)

plotMarkerComb2DScatters <- function(input, assigned_MarkerCombs, cell_subset_params = list(cellTypes = NULL,cellTypes_plot = NULL), palette_comb = NULL, return_list = FALSE){

  # extract marker names from the all positive population (tail of levels)
  marker_names <- strsplit(sub("\\+$", "", tail(levels(assigned_MarkerCombs),1)), "\\+")[[1]]

  # Determine expression matrix depending on input type
  if (methods::is(input, "matrix")) {
    # Already aggregated matrix
    exprs_matrix <- input

  } else if (methods::is(input, "SingleCellExperiment")) {
    # SCE input: extract and transpose expression assay
    exprs_matrix <- t(SummarizedExperiment::assay(input, "exprs")[marker_names, ])

  } else if (methods::is(input, "GatingSet") || methods::is(input, "flowSet")) {
    # Prepare flowSet
    fs <- if (methods::is(input, "GatingSet")) {
      flowWorkspace::cytoset_to_flowSet(flowWorkspace::gs_pop_get_data(input,"root"))
    } else {
      input
    }

    # Aggregate flowSet
    counts <- as.numeric(unlist(flowCore::fsApply(fs, flowCore::keyword, "$TOT"), use.names = FALSE))
    message("Aggregating data")
    agg <- suppressMessages(suppressWarnings(
      FlowSOM::AggregateFlowFrames(
        fs,
        channels = marker_names,
        cTotal = max(counts) * length(counts),
        writeOutput = FALSE,
        keepOrder = TRUE,
        sampleWithReplacement = FALSE
      )
    ))
    message("Done!")

    # Extract marker and channel names once
    param_data <- agg@parameters@data
    param_filtered <- param_data %>%
      dplyr::filter(.data$desc %in% marker_names)

    exprs_matrix <- agg@exprs[, param_filtered$name, drop = FALSE] %>%
      magrittr::set_colnames(param_filtered$desc)

  } else {
    stop("Unsupported input type. Expected matrix of cytokine expression, SingleCellExperiment, GatingSet, or flowSet.")
  }

  # construct dataframe to plot
  df_plot <- data.frame(exprs_matrix,
                        Polyfunctionality = assigned_MarkerCombs)

  subset_result <- .subsetDataCellTypes(df_plot,cell_subset_params = cell_subset_params)
  df_plot <- subset_result$dataframe
  plot_caption <- cowplot::ggdraw() + cowplot::draw_label(subset_result$plot_caption,size=8)

  df_CombinationSizes <- as.data.frame(table(df_plot$Polyfunctionality))

  # make long format for plotting
  df_plot_long <- df_plot %>%
    tidyr::pivot_longer(cols = c(-.data$Polyfunctionality),
                        names_to = "marker",
                        values_to = "exprs")%>%
    dplyr::mutate(marker = factor(.data$marker,levels = marker_names),
                  marker_num = as.numeric(.data$marker))%>%
    dplyr::group_by(.data$Polyfunctionality,.data$marker)%>%
    dplyr::mutate(scaled_y = if (dplyr::n() > 1) scale(.data$exprs)[, 1] else 0,
                  base_jitter = stats::runif(dplyr::n(), -0.15, 0.15),
                  density_scaled_jitter = if (dplyr::n() > 1) stats::rnorm(dplyr::n(), mean = 0, sd = 0.25 * stats::dnorm(.data$scaled_y)) else 0,
                  x_blob = .data$marker_num[1] + .data$base_jitter + .data$density_scaled_jitter) %>%
    dplyr::ungroup()

  if(is.null(palette_comb)){
    palette_comb <- .getNamedPaletteMarkerComb(marker_names = marker_names)
  }

  list_plots <- list()
  for(comb_plot in names(palette_comb)){

    color_highlight <- palette_comb[[comb_plot]]
    plot_title <- .markerCombtoMarkdown(comb_plot)
    n_pos_cells <- df_CombinationSizes %>%
      dplyr::filter(.data$Var1==comb_plot) %>%
      dplyr::pull(.data$Freq)

    p <- df_plot_long %>%
      dplyr::mutate(Color = .data$Polyfunctionality == comb_plot,
                    Polyfunctionality = factor(.data$Polyfunctionality, levels = c(comb_plot,setdiff(names(palette_comb),comb_plot)))) %>%
      dplyr::group_by(.data$Polyfunctionality) %>%
      dplyr::group_modify(~ {
        if (.y$Polyfunctionality != comb_plot & nrow(.x)>1000) { # downsample all background groups, i.e. grey points
          dplyr::sample_n(.x, 1000, replace = FALSE)
        } else {
          .x
        }
      }) %>%
      dplyr::ungroup()%>%
      dplyr::arrange(dplyr::desc(.data$Polyfunctionality))%>%
      ggplot(aes(x = .data$x_blob,
                 y = .data$exprs,
                 color = .data$Color))+
      geom_point(size = .8)+
      xlab(NULL)+
      ylab(NULL)+
      scale_x_continuous(breaks = sort(unique(df_plot_long$marker_num)),
                         labels = sort(unique(df_plot_long$marker)),
                         expand = expansion(mult = c(0.05, 0.05))) +
      labs(title = plot_title,
           subtitle = paste("n =",n_pos_cells))+
      theme_classic()+
      theme(plot.title = ggtext::element_markdown(hjust = .5,
                                                  face = "bold"),
            plot.subtitle = element_text(hjust = .5,
                                         size = 10),
            legend.position = "none",
            panel.grid.major = element_line(linewidth = 0.08,
                                            linetype = 'solid',
                                            colour = "grey60"),
            axis.text.x = element_text(angle = 90,
                                       hjust = 1,
                                       vjust = .5,
                                       size = 12))+
      scale_color_manual(values = c("TRUE" = color_highlight,
                                    "FALSE" = "grey85"))


    list_plots[[as.character(comb_plot)]] <- p

  }

  if(return_list == TRUE){

    return(list_plots)

  }else{

    p_combined <- cowplot::plot_grid(plotlist = list_plots)
    p_final <- cowplot::plot_grid(p_combined,plot_caption,ncol=1, rel_heights = c(1,.05))
    return(p_final)

  }


}

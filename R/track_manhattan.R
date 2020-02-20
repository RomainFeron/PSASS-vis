#' @title Plot a manhattan track
#'
#' @description Plot a single track for a manhattan plot
#'
#' @param data Input data generated with the \code{\link{create_manhattan_track_data}} function)
#'
#' @param track Track object storing properties for the current track, generated with the \code{\link{manhattan_track}} function
#'
#' @param bottom.track If TRUE, x-axis labels and title will be added to the plot (default: TRUE)
#'
#' @param show.chromosome.names If TRUE, display chromosome names on the x axis (default: TRUE)
#'
#' @param x.labels.angle If TRUE, replace chromosome names with numbers for readability (default: FALSE)
#'
#' @param x.axis.title Title to display on the x axis (default: NULL, i.e. no title)
#'
#' @return A ggplot object for the plot
#'
#' @examples
#' genomic_data <- load_genome_input("psass_window.tsv")
#' fst_track <- manhattan_track("Fst")
#' manhattan_data <- create_manhattan_track_data(genomic_data, fst_track)
#'
#' fst_plot <- plot_track_manhattan(manhattan_data, fst_track, bottom.track=TRUE)
#'

plot_track_manhattan <- function(data, backgrounds, track, bottom.track = TRUE,
                                 show.chromosome.names = TRUE, x.labels.angle = 90, x.axis.title = NULL) {


    if (is.null(x.axis.title)) { x.axis.title <- ggplot2::element_blank()}

    # Maximum / minimum Y value
    if (is.null(track$ylim)) { track$ylim = c(min(data[, 3]), 1.025 * max(data[, 3]) + 0.01) }
    ymin <- track$ylim[1]
    ymax <- track$ylim[2]

    # Create a fake color palette for both backgrounds and points
    merged_color_palette <- setNames(c(unique(data$Color), unique(backgrounds$color)), c(unique(data$Color), unique(backgrounds$color)))

    manhattan_plot <- ggplot2::ggplot() +
        # Backgrounds with alternating colors
        ggplot2::geom_rect(data = backgrounds,
                           ggplot2::aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = color),
                           alpha = 0.5) +
        # Data points
        ggplot2::geom_point(data = data,
                            ggplot2::aes(x = Position_plot, y = Value, color = Color),
                            size = track$point.size,
                            alpha = 1) +
        # Attribute color values from merged color scale for points and backgrounds
        ggplot2::scale_color_manual(values = merged_color_palette) +
        ggplot2::scale_fill_manual(values = merged_color_palette) +
        # Generate y-axis
        ggplot2::scale_y_continuous(name = track$label,
                                    limits = c(ymin, ymax),
                                    expand = c(0, 0)) +
        # Adjust theme elements
        cowplot::theme_cowplot() +
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       axis.title.x = ggplot2::element_text(face = "bold", margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(face = "bold", margin = ggplot2::margin(0, 10, 0, 0)),
                       axis.text.y = ggplot2::element_text(color = "black", face = "bold"),
                       axis.text.x = ggplot2::element_text(color = "black", face = "bold", vjust = 0.5, angle = x.labels.angle))

    if (show.chromosome.names) {

        # Generate x-axis, use background start and end to place LG labels
        manhattan_plot <- manhattan_plot + ggplot2::scale_x_continuous(name = x.axis.title,
                                                                       breaks = backgrounds$start + (backgrounds$end - backgrounds$start) / 2,
                                                                       labels = unique(data$Contig_plot),
                                                                       expand = c(0, 0))

    } else {

        # Generate empty x-axis
        manhattan_plot <- manhattan_plot +
            ggplot2::scale_x_continuous(expand = c(0, 0), name = x.axis.title) +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

    }
}
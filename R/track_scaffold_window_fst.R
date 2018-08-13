#' @title Scaffold FST window track
#'
#' @description Draws a scaffold plot track for FST window data.
#' This function is intended for use in the \code{\link{draw_scaffold_plot}} function.
#'
#' @param data FST window data frame.
#'
#' @param scaffold.name Name of the plotted scaffold (for the x axis)
#'
#' @param region A vector specifying the boundaries of the region to be plotted, e.g. c(125000, 250000).
#'
#' @param major.lines.y If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param major.lines.x If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param ylim Limits of the y axis (default: c(min(FST), 1)).
#'
#' @param color Color of the plotted area (default: "grey50").
#'
#' @param bottom.track If TRUE, this track will be considered bottom track of the plot and the x axis will be drawn (default: FALSE).


track_scaffold_window_fst <- function(data,
                                      scaffold.name,
                                      region,
                                      major.lines.y = TRUE,
                                      major.lines.x = FALSE,
                                      ylim = c(0, max(data$Fst)),
                                      color = "grey50",
                                      bottom.track = FALSE) {

    # Create major grid lines for y axis if specified
    if (major.lines.y) {
        major_lines_y <- ggplot2::element_line(color = "grey90", linetype = 1)
    } else {
        major_lines_y <- ggplot2::element_blank()
    }

    # Create major grid lines for x axis if specified
    if (major.lines.x) {
        major_lines_x <- ggplot2::element_line(color = "grey90", linetype = 1)
    } else {
        major_lines_x <- ggplot2::element_blank()
    }

    # Add x axis if bottom track
    if (!bottom.track) {
        axis_title_x <- ggplot2::element_blank()
    } else {
        axis_title_x <- ggplot2::element_text()
    }

    # Draw the plot
    g <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        ggplot2::geom_ribbon(data = data, ggplot2::aes(x = Original_position, ymin = 0, ymax = Fst), fill = color, color = "black") +
        ggplot2::scale_y_continuous(name = expression(paste("F"["ST"], " win.", sep="")), expand = c(0.01, 0.01), limits = ylim) +
        generate_x_scale(region, scaffold.name) +
        ggplot2::theme(legend.position = "none",
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(l = 5)),
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x,
                       axis.title.x = axis_title_x)

    return(g)
}

#' @title Scaffold combined depth window track
#'
#' @description Draws a scaffold plot track for both sexes depth window data.
#' This function is intended for use in the \code{\link{draw_scaffold_plot}} function.
#'
#' @param data depth window data frame.
#'
#' @param scaffold.name Name of the plotted scaffold (for the x axis)
#'
#' @param region A vector specifying the boundaries of the region to be plotted, e.g. c(125000, 250000).
#'
#' @param type Type of depth to use, either "absolute" or "relative" (default: "absolute").
#'
#' @param major.lines.y If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param major.lines.x If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param ylim Limits of the y axis (default: NULL).
#'
#' @param males.color Color of the plotted area for males (default: "dodgerblue3").
#'
#' @param females.color Color of the plotted area for females (default: "firebrick2").
#'
#' @param bottom.track If TRUE, this track will be considered bottom track of the plot and the x axis will be drawn (default: FALSE).
#'
#' @param alpha.value Transparency values for the tracks (default: 0.25).


track_scaffold_window_depth_combined <- function(data, scaffold.name, region, type = "absolute",
                                                 major.lines.y = TRUE, major.lines.x = FALSE,
                                                 ylim = NULL, males.color = "dodgerblue3", females.color = "firebrick2",
                                                 bottom.track = FALSE, alpha.value = 0.25) {

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

    # Create data based on type of depth
    if (type == "absolute") {
        data$Males <- data$Males_depth_abs
        data$Females <- data$Females_depth_abs
    } else if (type == "relative") {
        data$Males <- data$Males_depth_rel
        data$Females <- data$Females_depth_rel
    } else {
        stop(paste0(" - Error: depth type \"", type, "\" does not exist."))
    }

    # Y axis limits
    ylim <- c(0, max(data$Males, data$Females))

    # Add x axis if bottom track
    if (!bottom.track) {
        axis_title_x <- ggplot2::element_blank()
    } else {
        axis_title_x <- ggplot2::element_text()
    }

    # Draw the plot
    g <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        ggplot2::geom_ribbon(data = data, ggplot2::aes(x = Original_position, ymin = 0, ymax = Males, fill = "Males", color = "Males"),
                             alpha = alpha.value, show.legend = c("color"=TRUE, "fill"=FALSE)) +
        ggplot2::geom_ribbon(data = data, ggplot2::aes(x = Original_position, ymin = 0, ymax = Females, fill = "Females", color = "Females"),
                             alpha = alpha.value, show.legend = c("color"=TRUE, "fill"=TRUE)) +
        ggplot2::scale_y_continuous(name = "Combined depth window", expand = c(0.01, 0.01), limits = ylim) +
        generate_x_scale(region, scaffold.name) +
        ggplot2::scale_fill_manual(name = "", values = c("Males" = males.color, "Females" = females.color)) +
        ggplot2::scale_color_manual(name = "", values = c("Males" = males.color, "Females" = females.color)) +
        ggplot2::theme(legend.position = "top",
                       legend.margin = ggplot2::margin(t=10, r=0, b=0, l=10),
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(l = 5)),
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x,
                       axis.title.x = axis_title_x)
    return(g)
}

#' @title Scaffold depth window track
#'
#' @description Draws a scaffold plot track for single sex SNP depth data.
#' This function is intended for use in the \code{\link{draw_scaffold_plot}} function.
#'
#' @param data depth window data frame.
#'
#' @param sex Sex to plot the data for, either "males" or "females".
#'
#' @param scaffold.name Name of the plotted scaffold (for the x axis).
#'
#' @param type Type of depth to use, either "absolute" or "relative" (default: "absolute").
#'
#' @param region A vector specifying the boundaries of the region to be plotted, e.g. c(125000, 250000).
#'
#' @param major.lines.y If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param major.lines.x If TRUE, major grid lines will be plotted for the y axis (default: TRUE).
#'
#' @param ylim Limits of the y axis (default: NULL).
#'
#' @param color Color of the plotted area. If NULL, "dodgerblue3" will be used for males and "firebrick2" for females (default: NULL).
#'
#' @param bottom.track If TRUE, this track will be considered bottom track of the plot and the x axis will be drawn (default: FALSE).


track_scaffold_window_depth <- function(data, sex, scaffold.name, region, type = "absolute",
                                        major.lines.y = TRUE, major.lines.x = FALSE,
                                        ylim = NULL, color = NULL, bottom.track = FALSE) {

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

    # Generate data according to specified depth type
    if (type == "absolute") {
        data$Males <- data$Males_depth_abs
        data$Females <- data$Females_depth_abs
    } else if (type == "relative") {
        data$Males <- data$Males_depth_rel
        data$Females <- data$Females_depth_rel
    } else {
        stop(paste0(" - Error: depth type \"", type, "\" does not exist."))
    }

    # Generate data according to specified sex
    if (sex == "males") {
        cov_data <- data$Males
        if (is.null(color)) color <- "dodgerblue3"
        if (is.null(ylim)) ylim <- c(0, max(data$Males))
        sex_short <- "M."
    } else if (sex == "females") {
        cov_data <- data$Females
        if (is.null(color)) color <- "firebrick2"
        if (is.null(ylim)) ylim <- c(0, max(data$Females))
        sex_short <- "F."
    } else {
        stop(paste0(" - Error: sex \"", sex, "\" does not exist."))
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
        ggplot2::geom_ribbon(data = data, ggplot2::aes(x = Original_position, ymin = 0, ymax = cov_data),
                             fill = color, color = color, size = 0.4, alpha = 0.75) +
        ggplot2::scale_y_continuous(name = paste0(sex_short, " depth window"), expand = c(0.01, 0.01), limits = ylim) +
        generate_x_scale(region, scaffold.name) +
        ggplot2::theme(legend.position = "none",
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(l = 5)),
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x,
                       axis.title.x = axis_title_x)

    return(g)
}

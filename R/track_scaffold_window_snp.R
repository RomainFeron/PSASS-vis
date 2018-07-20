#' @title Scaffold SNP window track
#'
#' @description Draws a scaffold plot track for single sex SNP window data.
#' This function is intended for use in the \code{\link{draw_scaffold_plot}} function.
#'
#' @param data SNP window data frame.
#'
#' @param sex Sex to plot the data for, either "males" or "females".
#'
#' @param scaffold.name Name of the plotted scaffold (for the x axis)
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


track_scaffold_window_snp <- function(data,
                                      sex,
                                      scaffold.name,
                                      region,
                                      major.lines.y = TRUE,
                                      major.lines.x = FALSE,
                                      ylim = NULL,
                                      color = NULL,
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

    if (sex == "males") {
        snp_data <- data$Males
        if (is.null(color)) color <- "dodgerblue3"
        if (is.null(ylim)) ylim <- c(0, max(data$Males))
        sex_short <- "M."
    } else if (sex == "females") {
        snp_data <- data$Females
        if (is.null(color)) color <- "firebrick2"
        if (is.null(ylim)) ylim <- c(0, max(data$Females))
        sex_short <- "F."
    } else {
        stop(paste0(" - Error: sex \"", sex, "\" does not exist."))
    }

    # Draw the plot
    g <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        ggplot2::geom_ribbon(data = data, ggplot2::aes(x = Position, ymin = 0, ymax = snp_data),
                             fill = color, color = color, size = 0.4, alpha = 0.75) +
        ggplot2::scale_y_continuous(name = paste0(sex_short, " SNP window"), expand = c(0.01, 0.01), limits = ylim) +
        generate_x_scale(region, scaffold.name) +
        ggplot2::theme(legend.position = "none",
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(l = 5)),
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x)

    # Add x axis if bottom track
    if (!bottom.track) {
        g <- g + ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank())
    }

    return(g)
}

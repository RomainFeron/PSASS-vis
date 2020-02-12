track_region <- function(data, region_info, track, bottom.track = FALSE) {

    # Create major grid lines for y axis if specified
    if (track$major.lines.y) {
        major_lines_y <- ggplot2::element_line(color = "grey95", linetype = 1)
    } else {
        major_lines_y <- ggplot2::element_blank()
    }

    # Create major grid lines for x axis if specified
    if (track$major.lines.x) {
        major_lines_x <- ggplot2::element_line(color = "grey95", linetype = 1)
    } else {
        major_lines_x <- ggplot2::element_blank()
    }

    # Add x axis if bottom track
    if (!bottom.track) {
        axis_title_x <- ggplot2::element_blank()
    } else {
        axis_title_x <- ggplot2::element_text()
    }

    n_datasets <- (ncol(data) - 1) / 2

    # Check if at least one dataset is to be plotted as ribbon to setup y-axis limits
    has_ribbon = TRUE
    for (i in 1:n_datasets) { if (track$type[i] == "ribbon") { has_ribbon <- TRUE }}

    # Create y-axis limits if not specified, expand a bit from min and max
    if (is.null(track$ylim)) {

        ymin <- min(data[, 2 * seq(1, n_datasets)])
        if (has_ribbon == TRUE) ymin <- min(0, ymin)
        ymax <- max(data[, 2 * seq(1, n_datasets)])
        track$ylim <- c(0.975 * ymin, 1.025 * ymax)

    }

    # Draw the plot
    g <- ggplot2::ggplot() +
        cowplot::theme_cowplot() +
        ggplot2::scale_y_continuous(name = track$label, expand = ggplot2::expand_scale(c(0, 0.01), 0), limits = track$ylim) +
        generate_x_scale(region_info) +
        ggplot2::theme(legend.position = track$legend.position,
                       axis.text.y = ggplot2::element_text(margin = ggplot2::margin(l = 5)),
                       panel.grid.major.y = major_lines_y,
                       panel.grid.major.x = major_lines_x,
                       axis.title.x = axis_title_x)

    for (i in c(1:n_datasets)) {

        if (track$type[i] == "ribbon") {

            plot_data <- data[, c(1, 2*i)]
            names(plot_data) <- c("Position", "Values")
            ribbon_color = as.character(data[1, 2*i+1])

            g <- g + ggplot2::geom_ribbon(data = plot_data, ggplot2::aes(x = Position, ymin = 0, ymax = Values),
                                          fill = ribbon_color, color = ribbon_color, size = 0.4, alpha = track$alpha[i])

        } else if (track$type[i] == "points") {

            plot_data <- data[, c(1, 2*i, 2*i+1)]
            names(plot_data) <- c("Position", "Values", "Color")

            g <- g + ggplot2::geom_point(data = plot_data, ggplot2::aes(x = Position, y = Values),
                                         fill = Color, color = Color,
                                         size = point.size, alpha = track$alpha[i])

        }
    }

    return(g)
}
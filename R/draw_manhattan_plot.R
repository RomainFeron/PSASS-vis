#' @title Draw a manhattan plot
#'
#' @description Draw a manhattan plot from the PSASS data
#'
#' @param data A PSASS data structure obtained with the \code{\link{load_data_files}} function.
#'
#' @param output.file Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).
#'
#' @param width Width of the output file if specified, in inches (default: 14).
#'
#' @param height Height of the output file if specified, in inches (default: 8).
#'
#' @param dpi Resolution of the output file if specified, in dpi (default: 300).
#'
#' @param track Track to be plotted. Possible values are "position_fst", "window_fst", "window_snp_males", "window_snp_females", "depth_males", "depth_females"
#' (default: "window_fst").
#'
#' @param point.size Size of a point in the plot (default 0.5)
#'
#' @param point.palette Color palette for the dots (default c("dodgerblue3", "darkgoldenrod2"))
#'
#' @param background.palette Color palette for the background (default c("grey85", "grey100"))
#'
#' @param depth.type Type of depth to be plotted, either "absolute" or "relative" (default: "absolute").
#'
#' @param min.depth Minimum depth to compute depth ratio.
#' The ratio for positions with depth lower than this value in either sex will be 1 (default: 10).
#'
#' @examples
#'
#' draw_manhattan_plot(data, track = "window_fst")


draw_manhattan_plot <- function(data,
                                output.file = NULL, width = 14, height = 8, dpi = 300,
                                track = "window_fst",
                                point.size = 0.5, point.palette = c("dodgerblue3", "darkgoldenrod2"), background.palette = c("grey85", "grey100"),
                                depth.type = "absolute", min.depth = 10) {

    # Compute cumulative lengths of contigs, which will be added to the position of each point depending on the contig
    cumulative_lengths <- c(0, cumsum(data$lengths$lg))
    names(cumulative_lengths) <- c(names(data$lengths$lg), "Unplaced")

    # Load data
    if (track == "window_fst") {

        manhattan_data <- data$window_fst[, c(1, 2, 3, 5, 6)]
        y_title = expression(paste("F"["ST"], " in a sliding window"))

    } else if (track == "position_fst") {

        manhattan_data <- data$position_fst[, c(1, 2, 3, 5, 6)]
        y_title = expression(paste("Per-base F"["ST"]))

    } else if (track == "window_snp_males") {

        manhattan_data <- data$window_snp[, c(1, 2, 3, 6, 7)]
        y_title = expression(paste("Male-specific SNPs in a sliding window"))

    } else if (track == "window_snp_females") {

        manhattan_data <- data$window_snp[, c(1, 2, 4, 6, 7)]
        y_title = expression(paste("Female-specific SNPs in a sliding window"))

    } else if (track == "depth_males") {

        if (depth.type == "absolute") {

            manhattan_data <- data$depth[, c(1, 2, 3, 8, 9)]
            y_title = expression(paste("Absolute male depth in a sliding window"))

        } else if (depth.type == "relative") {

            manhattan_data <- data$depth[, c(1, 2, 5, 8, 9)]
            y_title = expression(paste("Relative male depth in a sliding window"))

        }

    } else if (track == "depth_females") {

        if (depth.type == "absolute") {

            manhattan_data <- data$depth[, c(1, 2, 4, 8, 9)]
            y_title = expression(paste("Absolute female depth in a sliding window"))

        } else if (depth.type == "relative") {

            manhattan_data <- data$depth[, c(1, 2, 6, 8, 9)]
            y_title = expression(paste("Relative female depth in a sliding window"))

        }

    }

    names(manhattan_data) <- c("Contig", "Position", "Value", "Contig_id", "Original_position")

    # Adjust x-axis position for each point in the data based on cumulative lengths of contigs
    manhattan_data$Position = manhattan_data$Position + cumulative_lengths[manhattan_data$Contig]

    # Attribute alternating colors to each contig
    order <- seq(1, length(data$lengths$plot))
    names(order) <- names(data$lengths$plot)
    manhattan_data$Color <- order[manhattan_data$Contig] %% 2
    manhattan_data$Color <- as.factor(as.character(manhattan_data$Color))

    # Create the background data for alternating background color
    background = data.frame(start = cumulative_lengths, end = cumulative_lengths + data$lengths$plot)
    background$Color = rep_len(c("A", "B"), length.out = dim(background)[1])  # Alternating A/B

    # Create merged color palette for background and data points
    merged_color_palette <- c("0"=point.palette[1], "1"=point.palette[2], "B"=background.palette[1], "A"=background.palette[2])

    # Maximum / minimum Y value
    ymax = 1.1 * max(manhattan_data$Value)
    ymin = min(manhattan_data$Value)

    manhattan_plot <- ggplot2::ggplot() +
        # Backgrounds with alternating colors
        ggplot2::geom_rect(data = background,
                           ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = ymax, fill = Color),
                           alpha = 0.5) +
        # Data points
        ggplot2::geom_point(data = manhattan_data,
                            ggplot2::aes(x = Position, y = Value, color = Color),
                            size = point.size,
                            alpha = 1) +
        # Attribute color values from merged color scale for points and backgrounds
        ggplot2::scale_color_manual(values = merged_color_palette) +
        ggplot2::scale_fill_manual(values = merged_color_palette) +
        # Generate x-axis, use background start and end to place LG labels
        ggplot2::scale_x_continuous(name="Linkage Group",
                                    breaks = background$start + (background$end - background$start) / 2,
                                    labels = names(order),
                                    expand = c(0, 0)) +
        # Generate y-axis
        ggplot2::scale_y_continuous(name=y_title,
                                    limits=c(ymin, ymax),
                                    expand = c(0, 0)) +
        # Adjust theme elements
        cowplot::theme_cowplot() +
        ggplot2::theme(legend.position = "none",
                       panel.grid.minor = ggplot2::element_blank(),
                       panel.grid.major.x = ggplot2::element_blank(),
                       axis.line.x = ggplot2::element_blank(),
                       panel.border = ggplot2::element_blank(),
                       axis.line.y = ggplot2::element_line(color = "black"),
                       axis.title.x = ggplot2::element_text(face="bold", margin = ggplot2::margin(10, 0, 0, 0)),
                       axis.title.y = ggplot2::element_text(face="bold", margin = ggplot2::margin(0, 10, 0, 0)),
                       axis.text.y = ggplot2::element_text(color="black", face="bold"),
                       axis.text.x = ggplot2::element_text(color="black", face="bold", angle=90, vjust=0.5))

    if (!is.null(output.file)) {

        cowplot::ggsave(output.file, plot = manhattan_plot, width = width, height = height, dpi = dpi)

    } else {

        print(manhattan_plot)

    }

    return(manhattan_plot)
}

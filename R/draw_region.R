
draw_region <- function(data, contig_lengths, region,
                        tracks = NULL,
                        output.file = NULL, width = 12, height = 4, res = 300) {

    # Add original contig names to contig lengths so that user can specify both chromosome names or contig names in region
    contig_lengths <- c(contig_lengths, setNames(unique(data$Length), unique(data$Contig)))

    # Get contig, min, and max from the region string
    region_info <- parse_region(region, contig_lengths)

    # Initialize list of plots
    n_tracks <- length(tracks)
    plots <- rep(list(NULL), n_tracks)

    # Draw specified tracks
    bottom_track <- FALSE
    for (i in c(1:n_tracks)) {

        if (i == n_tracks) bottom_track <- TRUE

        print(tracks[[i]]$label)

        track_data <- create_region_track_data(data, region_info, tracks[[i]])

        print(tracks[[i]]$label)
        print(tracks[[i]]$alpha)

        plots[[i]] <- track_region(track_data, region_info, tracks[[i]], bottom.track = bottom_track)

        print(tracks[[i]]$label)

    }

    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    print(combined)

    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    return(combined)
}



create_region_track_data <- function(data, region_info, track) {

    # Extract region from data
    if (region_info[[1]] %in% unique(data$Contig)) {

        data = subset(data, data$Contig == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    } else if (region_info[[1]] %in% unique(data$Contig_plot)) {

        data = subset(data, data$Contig_plot == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    }

    # Extract required columns and create color columns
    track_data <- data[, c("Position")]
    for (i in 1:length(track$metrics)) { track_data = cbind(track_data, data[, c(track$metrics[i])], rep(track$color[i], nrow(data))) }

    return(track_data)
}



track <- function(metrics, label = NULL, color = "grey20", alpha = 1, type = "ribbon", point.size = 0.5,
                  major.lines.y = TRUE, major.lines.x = FALSE, legend.position = "right", ylim=NULL) {

    n_metrics <- length(metrics)

    if (is.null(label)) {

        if (n_metrics ==1) {

            label <- metrics  # Set a default label if not specified

        } else {

            stop("Label required for multi-metrics track")

        }
    }

    if (n_metrics > 1) {

        # For each option, assign value to each metric if single value was defined in multi-metrics track
        if (length(color) == 1) { color <- rep(color, n_metrics) }
        if (length(alpha) == 1) { alpha <- rep(alpha, n_metrics) }
        if (length(type) == 1) { type <- rep(type, n_metrics) }
        if (length(point.size) == 1) { point.size <- rep(point.size, n_metrics) }

    }

    track_info <- list(metrics=metrics, label=label, color=color, alpha=alpha, type=type, point.size=point.size,
                       major.lines.y = major.lines.y, major.lines.x = major.lines.x, legend.position = legend.position,
                       ylim = ylim)

    return(track_info)
}

draw_region <- function(data, contig_lengths, region,
                        tracks = list("Fst", c("Snps_females", "Snps_males")),
                        track.labels = c("Fst", "Pool-specific SNPs"),
                        track.colors = list("grey20", c("firebrick2", "dodgerblue3")),
                        track.alpha = list(1, c(0.5, 0.5)),
                        track.types = list("ribbon", c("ribbon", "ribbon")),
                        point.size = 0.5, major.lines.y = TRUE, major.lines.x = FALSE,
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

        track_data <- create_region_track_data(data, region_info, tracks[[i]], track.colors[[i]])

        plots[[i]] <- track_region(track_data, region_info,
                                   plot.type = track.types[[i]], alpha = track.alpha[[i]], point.size = 0.1,
                                   major.lines.y = TRUE, major.lines.x = FALSE,
                                   ylim = NULL, bottom.track = bottom_track,
                                   y.label = track.labels[i], legend.position = "right")

    }

    combined <- cowplot::plot_grid(plotlist = plots, ncol = 1, align = "v")

    if (!is.null(output.file)) {

        ggplot2::ggsave(output.file, plot = combined, width = width, height = height * n_tracks, dpi = res)

    } else {

        print(combined)

    }

    return(combined)
}



create_region_track_data <- function(data, region_info, track, track.color) {

    if (region_info[[1]] %in% unique(data$Contig)) {

        data = subset(data, data$Contig == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    } else if (region_info[[1]] %in% unique(data$Contig_plot)) {

        data = subset(data, data$Contig_plot == region_info[[1]] & data$Position >= region_info[[2]] & data$Position <= region_info[[3]])

    }

    track_data <- data[, c("Position")]
    for (i in 1:length(track)) { track_data = cbind(track_data, data[, c(track[i])], rep(track.color[i], nrow(data))) }

    return(track_data)
}

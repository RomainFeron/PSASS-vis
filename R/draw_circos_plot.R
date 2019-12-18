#' @title Draw circos plot
#'
#' @description Draw a circos plot from the PSASS data
#'
#' @param data A PSASS data structure obtained with the \code{\link{load_data_files}} function.
#'
#' @param output.file Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).
#'
#' @param width Width of the output file if specified, in pixels (default: 2400).
#'
#' @param height Height of the output file if specified, in pixels (default: 2400).
#'
#' @param res Resolution of the output file if specified, in \% (default: 120).
#'
#' @param tracks Tracks to be plotted. Possible values are "fst_pos", "fst_win", "snp_pos", "snp_win_pool1",
#' "snp_win_pool2", "snp_win", "depth_pool1", "depth_pool2", "depth_ratio"
#' (default: c("fst_win", "snp_win_pool1", "snp_win_pool2", "depth_ratio")).
#'
#' @param highlight A vector of sectors to highlight, for instance c("LG5") or c("NC_02536.1", "NC_02543.1") (default: NULL).
#'
#' @param zoom.highlights If TRUE, highlighted sectors will be enlarged and placed at the top of the plot (default: FALSE).
#'
#' @param zoom.ratio Zoom factor for highlighted sectors if zoom.highlights is TRUE (default: 2).
#'
#' @param zoom.suffix Suffix to append to the name of zoomed highlighted sectors (default: " (zoom)").
#'
#' @param base.color Background color of a non-highlighted sector (default: "white").
#'
#' @param highlight.color Background color of a highlighted sector (default: "grey80").
#'
#' @param point.size Size of the points in the plot (default: 0.1).
#'
#' @param points.color Color of points for fst_pos, fst_win, and depth_ratio tracks (default: "grey20)
#'
#' @param pools.color Color of points for each pool for snp_pos, snp_win_pool1, snp_win_pool2, depth_pool1, and depth_pool2 tracks
#' (default: c("firebrick2", "dodgerblue3") for pool1 and pool2 respectively)
#'
#' @param color.unplaced If TRUE, unplaced scaffolds will be colored with alternating colors, like in a manhattan plot (default: FALSE)
#'
#' @param points.color.unplaced Alternating point color for scaffolds when plotting unplaced scaffolds (default: c(dodgerblue3", "goldenrod1"))
#'
#' @param sector.titles.expand Parameter to manually override the space between sector titles and x-axis (default: NULL).
#'
#' @param depth.type Type of depth to be plotted, either "absolute" or "relative" (default: "absolute").
#'
#' @param min.depth Minimum depth to compute depth ratio.
#' The ratio for positions with depth lower than this value in either sex will be 1 (default: 10).
#'
#' @examples
#'
#' draw_circos_plot(data, tracks = c("position_fst", "position_snp", "depth_males"),
#'                  highlight = c("Chr10", "Chr17"), color.unplaced = TRUE)


draw_circos_plot <- function(data,
                             output.file = NULL, width = 2400, height = 2400, res = 120,
                             tracks = c("fst_win", "snp_win_pool1", "snp_win_pool2", "depth_ratio"),
                             highlight = NULL, zoom.highlights = FALSE, zoom.ratio = 2, zoom.suffix = " (zoom)",
                             base.color = "white", highlight.color = "grey80", point.size = 0.1,
                             points.color = "grey20",
                             pools.color = c("firebrick2", "dodgerblue3"),
                             color.unplaced = FALSE,
                             points.color.unplaced = c("dodgerblue3", "goldenrod1"),
                             sector.titles.expand = NULL, depth.type = "absolute", min.depth = 10) {


    # Open output file if specified
    if (!is.null(output.file)) {
        png(output.file, width=width, height=height,  res=res)
    }


    # Get sector information
    n_sectors <- length(data$lengths$plot)
    sectors <- names(data$lengths$plot)
    sector_width <- data$lengths$plot
    if (!is.null(data$names)) {
        sector_names <- data$names
    } else {
        sector_names <- sectors
        names(sector_names) <- sectors
    }


    # Test if lgs specified in "highlight" exist, otherwise remove them from the highlight vector and print a warning
    if (!is.null(highlight)) {
        to_remove <- c()
        for (i in 1:length(highlight)) {
            if (!(highlight[i] %in% sectors)) {
                if (!is.null(data$names) & highlight[i] %in% data$names) {
                    highlight[i] <- names(data$names[which(data$names == highlight[i])])
                } else {
                    warning(paste0("Could not find sector \"", highlight[i], "\" given by parameter \"highlight\"."))
                    to_remove <- c(to_remove, i)
                }
            }
        }
        if (!is.null(to_remove)) {highlight <- highlight[-to_remove]}
    }


    # Create the zoomed sector if specified
    if (zoom.highlights & !is.null(highlight)) {

        zoom_data <- data$data[data$data$Contig %in% highlight, , drop = FALSE]  # Extract zoomed sectors data
        zoom_data$Contig <- paste0(zoom_data$Contig, zoom.suffix)  # Edit zoomed sector names with suffix
        data$data <- rbind(data$data, zoom_data)  # Combine base sectors data and zoomed sectors data
        zoom_lengths <- data$lengths$plot[highlight]  # Extract zoomed sector lengths
        names(zoom_lengths) <- paste0(highlight, zoom.suffix)  # Eddit zoomed sector names with suffix
        data$lengths$plot <- c(data$lengths$plot, zoom_lengths)  # Combine base sectors lengths and zoomed sectors lengths
        sectors <- c(sectors, paste0(highlight, zoom.suffix))  # Update sectors list
        n_sectors <- n_sectors + length(highlight)  # Update sectors count
        sector_width <- c(sector_width, zoom_lengths * zoom.ratio)  # Update sector widths
        zoom_names <- paste0(sector_names[highlight], zoom.suffix)
        names(zoom_names) <- paste0(highlight, zoom.suffix)  # Edit zoomed names with suffix
        sector_names <- c(sector_names, zoom_names)
    }


    # Setup sector colors
    bgs <- rep(base.color, n_sectors)
    bgs[which(sectors %in% highlight)] <- highlight.color
    bgs[which(sectors %in% paste0(highlight, zoom.suffix))] <- highlight.color


    # Setup point colors
    if (!color.unplaced) {
        point.color.unplaced = c(points.color, points.color)
    }


    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, n_sectors), data$lengths$plot), ncol=2)


    # Setup gaps between sectors
    gaps <- rep(1, n_sectors)
    gaps[length(gaps)] <- 12  # Bigger gap after the last sector to leave some space for track titles
    if (zoom.highlights & (length(highlight) > 0)) { gaps[length(gaps) - length(highlight)] <- 5 }  # Add a small gap before the first zoomed sector


    # Calculate angle offset to have zoomed sector on top
    if (zoom.highlights & length(highlight) > 0) {
        a <- 360 * sum(tail(sector_width, n = length(highlight))) / sum(sector_width) / 2 + 7.5
    } else {
        a <- 7.5
    }

    # Calculate width of track and sector.title.expand factor based on number of tracks
    n_tracks <- length(tracks)
    track.height <- 0
    if (is.null(sector.titles.expand)) sector.titles.expand <- 1.1 + n_tracks / 10
    if (n_tracks == 1) {
        track.height <- 0.3
    } else if (n_tracks >= 2) {
        track.height <- 0.5 / n_tracks
    }


    # Reset circos parameters
    circlize::circos.clear()


    # Setup circos parameters
    circlize::circos.par("track.height" = track.height,
                         "start.degree" = 90 - a,
                         "gap.degree" = gaps,
                         "cell.padding" = c(0.01, 0.1, 0.01, 0.1),
                         "points.overflow.warning" = FALSE,
                         "track.margin" = c(0, 0.015),
                         "unit.circle.segments" = 100)


    # Initialize circos plot
    circlize::circos.initialize(factors = sectors, xlim = sector_lengths, sector.width = sector_width)


    # Draw specified tracks
    for (i in c(1:length(tracks))) {

        top.track <- FALSE
        if (i == 1) top.track <- TRUE

        if (tracks[i] == "fst_pos") {

            track_circos_position_fst(data$position_fst, bg.col = bgs,
                                      point.size = point.size, top.track = top.track,
                                      points.color = points.color, points.color.unplaced = points.color.unplaced,
                                      sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "fst_win") {

            track_circos_window_fst(data$window_fst, bg.col = bgs,
                                    point.size = point.size, top.track = top.track,
                                    points.color = points.color, points.color.unplaced = points.color.unplaced,
                                    sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "snp_pos") {

            # TO BE IMPLEMENTED

            print("snp_pos track is not implemented yet")

        } else if (tracks[i] == "snp_win_pool1") {

            track_circos_window_snp(data$window_snp, pool = 1,
                                    bg.col = bgs, point.size = point.size, top.track = top.track,
                                    sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                                    pools.color = pools.color)

        } else if (tracks[i] == "snp_win_pool2") {

            track_circos_window_snp(data$window_snp, pool = 2,
                                    bg.col = bgs, point.size = point.size, top.track = top.track,
                                    sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                                    pools.color = pools.color)

        }  else if (tracks[i] == "snp_win") {

            # TO BE IMPLEMENTED

            print("snp_win track is not implemented yet")

            # track_circos_window_snp(data$window_snp, sex = "Females", bg.col = bgs,
            #                 top.track = top.track, sector.names = sector_names,
            #                 sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "depth_pool1") {

            track_circos_depth(data$depth, pool = 1, type = depth.type,
                               bg.col = bgs, top.track = top.track, point.size = point.size,
                               sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                               pools.color = pools.color)

        } else if (tracks[i] == "depth_pool2") {

            track_circos_depth(data$depth, pool = 2, type = depth.type,
                               bg.col = bgs, top.track = top.track, point.size = point.size,
                               sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                               pools.color = pools.color)

        } else if (tracks[i] == "depth_ratio") {

            track_circos_depth_ratio(data$depth, min.depth = min.depth,
                                     bg.col = bgs, point.size = point.size, top.track = top.track,
                                     points.color = points.color, points.color.unplaced = points.color.unplaced,
                                     sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else {

            print(paste0("Warning: unknown track type \"", tracks[i], "\" ..."))

        }
    }

    # Close output file
    if (!is.null(output.file)) {
        dev.off()
    }
}

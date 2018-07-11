#' @title Base circos plot
#'
#' @description Draw a circos plot from the PoolSex data
#'
#' @param data A PoolSex data structure obtained with the \code{\link{load_data_files}} function.
#'
#' @param tracks Tracks to be plotted. Possible values are "position_fst", "window_fst", "position_snp", "window_snp_males",
#' "window_snp_females", "combined_snp", "coverage_males", "coverage_females", "coverage_ratio"
#' (default c("window_fst", "combined_snp", "coverage_ratio")).
#'
#' @param highlight A vector of sectors to highlight, for instance c("LG5") or c("NC_02536.1", "NC_02543.1") (default NULL).
#'
#' @param zoom.highlights If TRUE, highlighted sectors will be enlarged and placed at the top of the plot (default FALSE).
#'
#' @param zoom.ratio Zoom factor for highlighted sectors if zoom.highlights is TRUE (default 2).
#'
#' @param zoom.suffix Suffix to append to the name of zoomed highlighted sectors (default " (zoom)").
#'
#' @param base.color Background color of a non-highlighted sector (default "white").
#'
#' @param highlight.color Background color of a highlighted sector (default "grey80").
#'
#' @param point.size Size of the points in the plot (default 0.1).
#'
#' @param color.unplaced If TRUE, unplaced scaffolds will be colored with alternating colors, like in a manhattan plot (default FALSE)
#'
#' @param color.palette Colors of the points in the plot. "0" and "1" specify the alternating colors for unplaced scaffolds
#' if color.unplaced is TRUE, "2" specifies the color for chromosomes for unsexed tracks (position_fst, window_fst, coverage_ratio), and
#' "males" and "females" specify the color of each sex for sexed tracks (position_snp, window_snp_males, window_snp_females, combined_snp,
#' coverage_males, coverage_females) (default c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20", "males"="dodgerblue3", "females"="firebrick2")).
#'
#' @param sector.titles.expand Parameter to control the space between sector titles and x-axis (default 1.5).
#'
#' @param coverage.type Type of coverage to be plotted, either "absolute" or "relative" (default "absolute").
#'
#' @param min.coverage Minimum coverage to compute coverage ratio.
#' The ratio for positions with coverage lower than this value in either sex will be 1 (default 10).
#'
#' @examples
#'
#' base_circos_plot(data, tracks = c("position_fst", "position_snp", "coverage_males"),
#'                  highlight = c("Chr10", "Chr17"), color.unplaced = TRUE)


base_circos_plot <- function(data,
                            tracks = c("window_fst", "window_snp", "coverage_ratio"),
                            highlight = NULL, zoom.highlights = FALSE, zoom.ratio = 2, zoom.suffix = " (zoom)",
                            base.color = "white", highlight.color = "grey80", point.size = 0.1,
                            color.unplaced = FALSE,
                            color.palette = c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20", "males"="dodgerblue3", "females"="firebrick2"),
                            sector.titles.expand = 1.5, coverage.type = "absolute", min.coverage = 10) {


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
        color.palette["0"] = color.palette["2"]
        color.palette["1"] = color.palette["2"]
    }


    # Create sector lengths matrix
    sector_lengths <- matrix(c(rep(0, n_sectors), data$lengths$plot), ncol=2)


    # Setup gaps between sectors
    gaps <- rep(1, n_sectors)
    gaps[length(gaps)] <- 12  # Bigger gap after the last sector to leave some space for track titles
    if (zoom.highlights) { gaps[length(gaps) - length(highlight)] <- 6 }  # Add a small gap before the first zoomed sector


    # Calculate angle offset to have zoomed sector on top
    a <- 360 * sum(tail(sector_width, n = length(highlight))) / sum(sector_width) / 2 + 6


    # Calculate width of track based on number of tracks
    n_tracks <- length(tracks)
    track.height <- 0
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

        if (tracks[i] == "position_fst") {

            draw_position_fst(data$position_fst, bg.col = bgs,
                              point.size = point.size, top.track = top.track, color.palette = color.palette,
                              sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "window_fst") {

            draw_window_fst(data$window_fst, bg.col = bgs,
                            point.size = point.size, top.track = top.track, color.palette = color.palette,
                            sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "position_snp") {

            # TO BE IMPLEMENTED

            draw_window_snp(data$window_snp, sex = "Males", bg.col = bgs,
                            top.track = top.track, sector.names = sector_names,
                            sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "window_snp_males") {

            draw_window_snp(data$window_snp, sex = "Males",
                            bg.col = bgs, point.size = point.size, top.track = top.track,
                            sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                            males.color = as.character(color.palette["males"]), females.color = as.character(color.palette["females"]))

        } else if (tracks[i] == "window_snp_females") {

            draw_window_snp(data$window_snp, sex = "Females",
                            bg.col = bgs, point.size = point.size, top.track = top.track,
                            sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                            males.color = as.character(color.palette["males"]), females.color = as.character(color.palette["females"]))

        }  else if (tracks[i] == "combined_snp") {

            # TO BE IMPLEMENTED

            draw_window_snp(data$window_snp, sex = "Females", bg.col = bgs,
                            top.track = top.track, sector.names = sector_names,
                            sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else if (tracks[i] == "coverage_males") {

            draw_coverage(data$coverage, sex = "Males", type = coverage.type,
                          bg.col = bgs, top.track = top.track, point.size = point.size,
                          sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                          males.color = as.character(color.palette["males"]), females.color = as.character(color.palette["females"]))

        } else if (tracks[i] == "coverage_females") {

            draw_coverage(data$coverage, sex = "Females", type = coverage.type,
                          bg.col = bgs, top.track = top.track, point.size = point.size,
                          sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors,
                          males.color = as.character(color.palette["males"]), females.color = as.character(color.palette["females"]))

        } else if (tracks[i] == "coverage_ratio") {

            draw_coverage_ratio(data$coverage, min.cov = min.coverage,
                                bg.col = bgs, point.size = point.size, top.track = top.track, color.palette = color.palette,
                                sector.names = sector_names, sector.titles.expand = sector.titles.expand, sectors = sectors)

        } else {

            print(paste0("Warning: unknown track type \"", tracks[i], "\" ..."))

        }
    }
}

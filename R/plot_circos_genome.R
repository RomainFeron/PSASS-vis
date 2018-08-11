#' @title Plot genome circos
#'
#' @description Generate a circos plot of the entire genome from the results of PoolSex analysis.
#'
#' @param contig_lengths_file_path Path to a contig lengths file.
#'
#' @param chromosomes_names_file_path Path to a contig names file (default NULL).
#'
#' @param prefix Prefix (including full path) to a complete dataset. If prefix is specified, it will be override individual file specifications
#' such as "window_fst_file_path" (default NULL).
#'
#' @param window_fst_file_path Path to a FST window output file (default NULL).
#'
#' @param position_fst_file_path Path to a FST positions output file (default NULL).
#'
#' @param window_snps_file_path Path to a SNPs window output file (default NULL).
#'
#' @param position_snps_file_path Path to a SNPs positions output file (default NULL).
#'
#' @param coverage_file_path Path to a coverage output file (default NULL).
#'
#' @param plot.unplaced If TRUE, unplaced contigs will be plotted as a supercontig (default TRUE).
#'
#' @param output.file Path to an output file in PNG format. If NULL, the plot will be drawn in the default graphic device (default: NULL).
#'
#' @param width Width of the output file if specified, in pixels (default: 2400).
#'
#' @param height Height of the output file if specified, in pixels (default: 2400).
#'
#' @param res Resolution of the output file if specified, in \% (default: 120).
#'
#' @param tracks Tracks to be plotted. Possible values are "position_fst", "window_fst", "position_snp", "window_snp_males",
#' "window_snp_females", "combined_snp", "coverage_males", "coverage_females", "coverage_ratio"
#' (default: c("window_fst", "combined_snp", "coverage_ratio")).
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
#' @param color.unplaced If TRUE, unplaced scaffolds will be colored with alternating colors, like in a manhattan plot (default: FALSE)
#'
#' @param color.palette Colors of the points in the plot. "0" and "1" specify the alternating colors for unplaced scaffolds
#' if color.unplaced is TRUE, "2" specifies the color for chromosomes for unsexed tracks (position_fst, window_fst, coverage_ratio), and
#' "males" and "females" specify the color of each sex for sexed tracks (position_snp, window_snp_males, window_snp_females, combined_snp,
#' coverage_males, coverage_females) (default: c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20", "males"="dodgerblue3", "females"="firebrick2")).
#'
#' @param sector.titles.expand Parameter to manually override the space between sector titles and x-axis (default: NULL).
#'
#' @param coverage.type Type of coverage to be plotted, either "absolute" or "relative" (default: "absolute").
#'
#' @param min.coverage Minimum coverage to compute coverage ratio.
#' The ratio for positions with coverage lower than this value in either sex will be 1 (default: 10).
#'
#' @examples
#'
#' # Standard plot with default options
#' plot_genome_circos(contig_lengths_file_path = 'data/contig_lengths.tsv',
#'                    chromosomes_names_file_path = 'data/chromosomes_names.tsv',
#'                    prefix = 'data/poolsex_analysis',
#'                    output.file = 'figures/poolsex_genome.png')
#'
#' # Plot FST positions and male SNP window with highlight on NC_02456.3
#' plot_genome_circos(contig_lengths_file_path = 'data/contig_lengths.tsv', prefix = 'data/poolsex_analysis',
#'                    tracks = c("position_fst", "window_snp_males"), highlight = c("NC_02456.3"))


plot_genome_circos <- function(contig_lengths_file_path,
                               plot.unplaced = TRUE, prefix = NULL,
                               window_fst_file_path = NULL, position_fst_file_path = NULL,
                               window_snps_file_path = NULL, position_snps_file_path = NULL,
                               coverage_file_path = NULL, chromosomes_names_file_path = NULL,
                               output.file = NULL, width = 2400, height = 2400, res = 120,
                               tracks = c("window_fst", "window_snp_males", "window_snp_females", "coverage_ratio"),
                               highlight = NULL, zoom.highlights = FALSE, zoom.ratio = 2, zoom.suffix = " (zoom)",
                               base.color = "white", highlight.color = "grey80", point.size = 0.1,
                               color.unplaced = FALSE,
                               color.palette = c("0"="dodgerblue3", "1"="goldenrod1", "2"="grey20", "males"="dodgerblue3", "females"="firebrick2"),
                               sector.titles.expand = NULL, coverage.type = "absolute", min.coverage = 10) {


    # Load the data
    data <- load_data_files(contig_lengths_file_path,
                            prefix = prefix,
                            window_fst_file_path = window_fst_file_path,
                            position_fst_file_path = position_fst_file_path,
                            window_snps_file_path = window_snps_file_path,
                            position_snps_file_path = position_snps_file_path,
                            coverage_file_path = coverage_file_path,
                            chromosomes_names_file_path = chromosomes_names_file_path,
                            plot.unplaced = plot.unplaced)

    # Draw the plot
    draw_circos_plot(data,
                     output.file = output.file,
                     width = width,
                     height = height,
                     res = res,
                     tracks = tracks,
                     highlight = highlight,
                     zoom.highlights = zoom.highlights,
                     zoom.ratio = zoom.ratio,
                     zoom.suffix = zoom.suffix,
                     base.color = base.color,
                     highlight.color = highlight.color,
                     point.size = point.size,
                     color.unplaced = color.unplaced,
                     color.palette = color.palette,
                     sector.titles.expand = sector.titles.expand,
                     coverage.type = coverage.type,
                     min.coverage = min.coverage)
}

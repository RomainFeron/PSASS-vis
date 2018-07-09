#' @title Load a single data file
#'
#' @description Loads all data files from poolsex analysis and other files needed for the plots.
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
#' @return A list with the following elements:
#' - names : chromosomes names (if specified)
#' - lengths : contig lengths
#' - window_fst : sliding window fst data (if specified)
#' - position_fst : position fst data (if specified)
#' - window_snp : sliding window snp data (if specified)
#' - position_snp : position snp data (if specified)
#' - coverage : coverage data (if specified)
#'
#' @examples
#' data <- load_data_files(prefix = "data/poolseq_analysis", chromosomes_names_file_path = "data/chromosomes_names.tsv",
#'                         contig_lengths_file_path = "data/contig_lengths.tsv")
#'
#' data <- load_data_files(window_fst_file_path = "data/poolseq_analysis_window_fst.tsv",
#'                         window_snps_file_path = "data/poolseq_analysis_window_snps.tsv",
#'                         coverage_file_path = "data/poolseq_analysis_coverage.tsv",
#'                         contig_lengths_file_path = "data/contig_lengths.tsv",
#'                         plot.unplaced = FALSE)


load_data_files <- function(contig_lengths_file_path,
                            prefix = NULL,
                            window_fst_file_path = NULL, position_fst_file_path = NULL,
                            window_snps_file_path = NULL, position_snps_file_path = NULL,
                            coverage_file_path = NULL, chromosomes_names_file_path = NULL,
                            plot.unplaced = TRUE) {

    output <- list()

    print(" - Loading chromosomes names file")
    output$names <- load_chromosomes_names(chromosomes_names_file_path)

    print(" - Loading contig lengths file")
    output$lengths <- load_contig_lengths(contig_lengths_file_path, chromosomes_names = output$names)

    if (!is.null(prefix)) {

        window_fst_file_path <- paste0(prefix, "_window_fst.tsv")
        position_fst_file_path <- paste0(prefix, "_position_fst.tsv")
        window_snp_file_path <- paste0(prefix, "_window_snp.tsv")
        position_snp_file_path <- paste0(prefix, "_position_snp.tsv")
        coverage_file_path <- paste0(prefix, "_coverage.tsv")

    }

    if (!is.null(position_fst_file_path)) {

        print(" - Loading positions FST file")
        output$position_fst <- load_single_data_file(position_fst_file_path, output$lengths, plot.unplaced = plot.unplaced)

    }

    if (!is.null(window_fst_file_path)) {

        print(" - Loading sliding window FST file")
        output$window_fst <- load_single_data_file(window_fst_file_path, output$lengths, plot.unplaced = plot.unplaced)

    }

    if (!is.null(position_snp_file_path)) {

        print(" - Loading positions SNP file")
        output$position_snp <- load_single_data_file(position_snp_file_path, output$lengths, plot.unplaced = plot.unplaced)

    }

    if (!is.null(window_snp_file_path)) {

        print(" - Loading sliding window SNP file")
        output$window_snp <- load_single_data_file(window_snp_file_path, output$lengths, plot.unplaced = plot.unplaced)

    }

    if (!is.null(coverage_file_path)) {

        print(" - Loading coverage file")
        output$coverage <- load_single_data_file(coverage_file_path, output$lengths, plot.unplaced = plot.unplaced)

    }

    return(output)
}
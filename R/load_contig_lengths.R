#' @title Load a contig lengths file
#'
#' @description Loads a table of contig lengths. The function separates the lengths of chromosomes and unplaced contigs
#' If chromosomes names are provided, they will be used to determine the chromosomes. Otherwise, the function will attempt to detect
#' chromosomes based on contig names (chromosomes start with LG / Chr / NC).
#'
#' @param input_file_path Path to a contig lengths file.
#'
#' @param chromosomes A vector of chromosomes names obtained with the \code{\link{load_chromosomes_names}} function (default NULL).
#'
#' @return A list with the following elements:
#' - lg : lengths of contigs determined to be chromosomes
#' - unplaced : lengths of contigs determined to be unplaced
#' - plot : lengths of sectors to be plotted
#'
#' @examples
#' lengths <- load_contig_lengths("contig_lengths.tsv", chromosomes_names = chromosomes_names)



load_contig_lengths <- function(input_file_path, chromosomes_names = NULL, plot.unplaced = TRUE) {

    raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
    data <- raw_data$X2
    names(data) <- raw_data$X1

    output <- list()

    if (!is.null(chromosomes_names)) {  # If a chromosomes names file was provided, it is used to determine chromosomes

        output$lg <- subset(data, names(data) %in% names(chromosomes_names))
        output$lg <- output$lg[gtools::mixedorder(chromosomes_names[names(output$lg)])]
        output$unplaced <- subset(data, !(names(data) %in% names(chromosomes_names)))

    } else {  # Otherwise, try to determine chromosomes automatically (condition: start with LG / Ch / NC)

        output$lg <- subset(data, substr(names(data), 1, 2) %in% c("LG", "lg", "Lg", "Ch", "ch", "CH", "NC"))

        if (is.vector(output$lg) && length(output$lg) > 1) {

            # Usually mitochondria is also called NC_xxx. If one chromosome is > 100 times smaller than the average of all ohter chromosomes,
            # it is considered to be the mitochondria and is removed
            putative_mt <- min(output$lg)
            if (100 * putative_mt < mean(output$lg)) output$lg <- output$lg[output$lg != putative_mt]

            output$lg <- output$lg[gtools::mixedorder(names(output$lg))]  # Order chromosomes based on their ID

        } else if (length(output$lg) == 1) {
            output$lg <- output$lg[-1,]  # If there is only one LG, it is assumed it is a mitochondria
        }

        output$unplaced <- subset(data, !(substr(names(data), 1, 2) %in% c("LG", "lg", "Ch", "ch", "MT", "mt", "NC")))

    }

    output$unplaced <- sort(output$unplaced, decreasing = TRUE)  # Sort unplaced scaffolds by size

    if (plot.unplaced) {

        output$plot <- c(output$lg, "Unplaced" = sum(output$unplaced))

    } else {

        output$plot <- output$lg

    }

    return(output)
}



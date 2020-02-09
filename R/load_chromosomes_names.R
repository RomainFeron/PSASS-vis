#' @title Load chromosome names
#'
#' @description Loads chromosome names from a tabulated file
#'
#' @param input_file_path Path to the chromosome names file
#'
#' @return A vector with contigs as names and chromosome names as values, or NULL if there is no input file
#'
#' @examples
#' chromosomes <- load_chromosome_names("chromosomes_names.tsv")


load_chromosome_names <- function(input_file_path) {

    if (!is.null(input_file_path)) {

        raw_data <- suppressMessages(readr::read_delim(input_file_path, "\t", escape_double = FALSE, col_names = FALSE,  trim_ws = TRUE))
        data <- raw_data$X2
        names(data) <- raw_data$X1
        data <- gtools::mixedsort(data)  # Sort alphabetically
        return(data)

    } else {

        return(NULL)

    }
}


parse_region <- function(region, contig_lengths) {

    tmp = strsplit(region, ":")
    contig = tmp[[1]][1]

    if (!(contig %in% names(contig_lengths))) {
        print(paste0("Error: invalid contig in region <", region, ">"))
        exit(1)
    }

    if (length(tmp[[1]]) == 2) {

        tmp = strsplit(tmp[[1]][2], "-")
        start = as.numeric(tmp[[1]][1])
        end = as.numeric(tmp[[1]][2])

    } else {

        start = 0
        end = unname(contig_lengths[contig])

    }

    return(list(contig, start, end))
}


convert_to_mb <- function(x, n = 0) {
    round(x / 10^6, n)
}


convert_to_kb <- function(x, n = 0) {
    round(x / 10^3, n)
}


# Generate a nice scale for any interval
generate_x_scale <- function(region_info) {

    # Size of the region
    S <- (region_info[[3]] - region_info[[2]]) / 10
    if (S <= 0) stop(paste0("Error: the size of the region has to be > 0 (value: ", S * 10, ")"))

    # Find the order of magnitude of the region's size
    N <- floor(log(S, 10))
    N10 <- 10 ^ N

    # Generate a scale of 10 round values within the region
    scale <- seq(N10 * floor(region_info[[2]] / N10), N10 * ceiling(region_info[[3]] / N10), N10 * round(S / N10, 1))

    # Adjust the labels based on the size of the values (megabp or kilop)
    if (region_info[[2]] < 10 ^ 6 & region_info[[3]] < 10 ^ 6) {
        scale_labels <- round(scale / 10 ^ 3, 3 - N)
        bp_unit <- "K"
    } else {
        scale_labels <- round(scale / 10 ^ 6, 6 - N)
        bp_unit <- "M"
    }

    output <- ggplot2::scale_x_continuous(name = paste0("Position on ", region_info[[1]], " (", bp_unit, "bp)"),
                                          expand = c(0.01, 0.01),
                                          breaks = scale,
                                          labels = scale_labels,
                                          limits = c(min(scale[1], region_info[[2]]), max(tail(scale, 1), region_info[[3]])))

    return(output)
}
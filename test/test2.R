data = psass_window_chr$data
contig.lengths = psass_window_chr$lengths
tracks = list(manhattan_track("Fst"))
output.file = NULL
width = 14
track.height = 6
res = 300
chromosomes.as.numbers = FALSE
default.point.color = c("dodgerblue3", "darkgoldenrod2")
default.bg.color = c("grey85", "white")
default.point.size = 0.5
default.ylim = NULL


data <- track_data
backgrounds <- track_background_data
track <- tracks[[1]]
bottom.track <- TRUE

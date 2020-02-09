setwd("/home/romain/work/code/PSASS-vis/tests/")

###################################
########## DATA LOADING ###########
###################################

# Load chromosome names from file
chromosomes = load_chromosome_names("chromosomes.tsv")

load_single_data_file
load_data_files


draw_circos_plot
plot_genome_circos

draw_manhattan_plot
plot_manhattan

draw_scaffold_plot
plot_scaffold


track_circos_depth
track_circos_depth_ratio
track_circos_position_fst
track_circos_window_fst
track_circos_window_snp
track_circos_window_snp_combined
track_scaffold_position_fst
track_scaffold_window_depth
track_scaffold_window_depth_combined
track_scaffold_window_depth_ratio
track_scaffold_window_fst
track_scaffold_window_snp
track_scaffold_window_snp_combined
track_scaffold_window_snp_ratio

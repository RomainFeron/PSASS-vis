setwd("/home/romain/work/code/PSASS-vis/test/")

###################################
########## DATA LOADING ###########
###################################

# Load chromosome names from file
chromosomes = load_chromosome_names("chromosomes.tsv")

# Load psass window output
psass_window_chr = load_genome_input("psass_window.tsv", chromosomes)
psass_window_chr_detect = load_genome_input("psass_window.tsv", chromosomes=NULL, detect.chromosomes = TRUE)
psass_window_no_chr = load_genome_input("psass_window.tsv", detect.chromosomes = FALSE)
psass_window_no_unplaced = load_genome_input("psass_window_no_unplaced.tsv")

# Load psass snp output
psass_snp_chr = load_genome_input("psass_snps.tsv", chromosomes)
psass_snp_chr_detect = load_genome_input("psass_snps.tsv", chromosomes=NULL, detect.chromosomes = TRUE)
psass_snp_no_chr = load_genome_input("psass_snps.tsv", detect.chromosomes = FALSE)
psass_snp_no_unplaced = load_genome_input("psass_snps_no_unplaced.tsv")

# Load psass fst output
psass_fst_chr = load_genome_input("psass_fst.tsv", chromosomes)
psass_fst_chr_detect = load_genome_input("psass_fst.tsv", chromosomes=NULL, detect.chromosomes = TRUE)
psass_fst_no_chr = load_genome_input("psass_fst.tsv", detect.chromosomes = FALSE)
psass_fst_no_unplaced = load_genome_input("psass_fst_no_unplaced.tsv")

# Circos plot
draw_circos_plot(data, contig_lengths, output.file = "circos.png")

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

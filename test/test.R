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
draw_circos_plot(psass_window_chr$data, psass_window_chr$lengths, output.file = "circos.png")

# Region plot
region = draw_region(psass_window_chr$data, psass_window_chr$lengths, "Chr24",
                     tracks = list(track("Fst", label = expression("F"["ST"])),
                                   track(c("Snps_females", "Snps_males"), label = "Pool-specific SNPs", color = c("firebrick2", "dodgerblue3"), alpha=0.6),
                                   track(c("Abs_depth_females", "Abs_depth_males"), label = "Absolute depth", color = c("firebrick2", "dodgerblue3"), alpha=c(0.4, 0.4))),
                     output.file = "region.png", width = 12, track.height = 4, res = 300)


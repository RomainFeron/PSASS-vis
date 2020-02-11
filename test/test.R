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
draw_circos_plot(psass_window_chr$data, psass_window_chr$contig_lengths, output.file = "circos.png")

region = draw_region(psass_window_chr$data, psass_window_chr$lengths, "Chr01:0-5000000",
                     tracks = list("Fst", c("Snps_females", "Snps_males"), c("Rel_depth_females", "Rel_depth_males")),
                     track.labels = c("Fst", "Pool-specific SNPs", "Depth"),
                     track.colors = list("grey20", c("firebrick2", "dodgerblue3"), c("firebrick2", "dodgerblue3")),
                     track.alpha = list(1, c(0.5, 0.5), c(0.5, 0.5)),
                     track.types = list("ribbon", c("ribbon", "ribbon"), c("ribbon", "ribbon")),
                     point.size = 0.5, major.lines.y = TRUE, major.lines.x = FALSE,
                     output.file = "region.png", width = 12, height = 4, res = 300)


# Parse command-line arguments
args = commandArgs(trailingOnly=T)
chip_bam_file = args[1]
input_bam_file = args[2]
bin_dir = args[3]
z_thresh = as.numeric(args[4])

# Load packages and set library path environment variable
#library("Rsamtools", lib.loc=paste0(bin_dir, "/util/")) # ML zainstalowalem te paczki globlnie,
# korzystaniem z tego co jest w tych katalogach byloby baaardzo ryzykowane (zbudowane dla archaicznego R 3.2 bodajze ...)
#library("fastcluster", lib.loc=paste0(bin_dir, "/util/")) # ML j.w
#library("spp", lib.loc=paste0(bin_dir, "/util/")) # ML j.w. 
#.libPaths(paste0(bin_dir, "/util/")) # hmm nie wiem po co to jest ...

# Load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos='https://ftp.gwdg.de/pub/misc/cran/')
BiocManager::install("Rsamtools", ask=FALSE)
install.packages('fastcluster', repos='https://ftp.gwdg.de/pub/misc/cran/') 
install.packages("spp", repos='https://ftp.gwdg.de/pub/misc/cran/', dependencies=TRUE)


library("Rsamtools")
library("fastcluster")
library("spp")

# Read the BAM files for the ChIP sample and the input-control sample
chip_reads = read.bam.tags(chip_bam_file)
input_reads = read.bam.tags(input_bam_file)


# Estimate of the binding-peak separation distance, 
# cross-correlation profile, and tag-quality bin-acceptance information
binding_char = get.binding.characteristics(
    chip_reads, srange=c(270, 1000), bin=5, accept.all.tags=T)


# Select tags based on the binding characteristics
chip_data = select.informative.tags(chip_reads, binding_char)
input_data = select.informative.tags(input_reads, binding_char)


# Remove or restrict singular positions with extremely high tag count
# relative to their neighborhood
chip_filt = remove.local.tag.anomalies(chip_data)
input_filt = remove.local.tag.anomalies(input_data)

# Get window size
window_half_size = binding_char$whs

# Identify broad regions of enrichment for the specified scale 
# ML Drobna uwaga czasami program wyrzuca sie bo dostaje na dalszych chromosomach nie numeryczne wartosci do funkcji min/max, ewentualnie dostaje Inf/-Inf
# mozna wiec zamiast podawac info o wszystkich chromosomach (ponad 3000!) dac tylko te ktore nas interesuja 
# broad_clusters = get.broad.enrichment.clusters(chip_filt[1:25], input_filt[1:25], window.size=window_half_size, z.thr=z_thresh, ...)

broad_clusters = get.broad.enrichment.clusters(
    chip_filt, input_filt, window.size=window_half_size, z.thr=z_thresh,
    tag.shift=round(binding_char$peak$x / 2))

# Write output file of broad clusters
suffix = paste0(".spp.z", z_thresh, ".broadPeak")
out_file = gsub(".bam", suffix, chip_bam_file)
write.broadpeak.info(broad_clusters, out_file)

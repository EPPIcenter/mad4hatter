library(argparse)
library(dada2)

# Argument Parsing and Setup
parser <- ArgumentParser(description='Perform Dereplication')
parser$add_argument('--filtFs-path', type="character", required=TRUE, help="Path to filtered F reads", nargs="+")
parser$add_argument('--filtRs-path', type="character", required=TRUE, help="Path to filtered R reads", nargs="+")
parser$add_argument('--filter-metadata', type="character", required=TRUE, help="Path to filter metadata")
parser$add_argument('--dout', type="character", required=TRUE, help="Output directory")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--verbose', action="store_true", help="Verbose")

args <- parser$parse_args()
print(args)

# Load filter metadata
filter_metadata <- readRDS(args$filter_metadata)

if (!dir.exists(args$dout)) {
    dir.create(args$dout, recursive=TRUE)
}

# Dereplicate and save for each pair of filtFs and filtRs
for(i in seq_along(args$filtFs_path)) {
    
    # Perform Dereplication for F reads
    derepF <- derepFastq(args$filtFs_path[i], verbose = args$verbose)
    # Construct unique filenames for saving F reads
    outFileF <- file.path(args$dout, paste0(basename(args$filtFs_path[i]), "_dereped.RDS"))
    saveRDS(derepF, file=outFileF)
    
    # Perform Dereplication for R reads
    derepR <- derepFastq(args$filtRs_path[i], verbose = args$verbose)
    # Construct unique filenames for saving R reads
    outFileR <- file.path(args$dout, paste0(basename(args$filtRs_path[i]), "_dereped.RDS"))
    saveRDS(derepR, file=outFileR)

}

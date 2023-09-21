library(argparse)
library(dada2)

# Argument Parsing and Setup
parser <- ArgumentParser(description='Learn Error Models')
parser$add_argument('--filt-1', nargs='+', type="character", required=TRUE, help="Paths to directories with filtered F reads or merged reads")
parser$add_argument('--filt-2', nargs='+', type="character", required=FALSE, help="Paths to directories with filtered R reads")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--maxConsist', type="numeric", default=10, help="Maximum number of iterations for model convergence")
parser$add_argument('--rand', action='store_true', help="Randomize")
parser$add_argument('--dout', type="character", required=TRUE, help="Output directory")

args <- parser$parse_args()
print(args)

if (!dir.exists(args$dout)) {
    dir.create(args$dout, recursive=TRUE)
}

# Learn errors for 1 or Merged using all loaded files
error_model <- learnErrors(args$filt_1, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand)

# If Reverse reads are not provided, save the error model as err_model.RDS
if (is.null(args$filt_2)) {
    saveRDS(error_model, file.path(args$dout, "err_model.RDS"))
} else {
    saveRDS(error_model, file.path(args$dout, "err_F_model.RDS"))
}

if (!is.null(args$filt_2)) {
    # Learn errors for 2 using all loaded files
    error_model <- learnErrors(args$filt_2, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand)
    saveRDS(error_model, file.path(args$dout, "err_R_model.RDS"))
}

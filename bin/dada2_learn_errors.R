library(logger)
load_library <- function(library_name) {
  output <- capture.output({
    suppressWarnings({
      library(library_name, character.only = TRUE)
    })
  }, type = "message")
  
  # Separate warnings from messages
  warnings <- warnings()
  
  # Log messages
  if(length(output) > 0) {
    log_info(paste("Message from", library_name, ":", paste(output, collapse = "; ")))
  }
  
  # Log warnings
  if(length(warnings) > 0) {
    log_warn(paste("Warning from", library_name, ":", paste(warnings, collapse = "; ")))
  }
}

load_library("argparse")
load_library("dada2")

# Argument Parsing and Setup
parser <- ArgumentParser(description='Learn Error Models')
parser$add_argument('--filt-1', nargs='+', type="character", required=TRUE, help="Paths to directories with filtered F reads or merged reads")
parser$add_argument('--filt-2', nargs='+', type="character", required=FALSE, help="Paths to directories with filtered R reads")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--maxConsist', type="numeric", default=10, help="Maximum number of iterations for model convergence")
parser$add_argument('--rand', action='store_true', help="Randomize")
parser$add_argument('--dout', type="character", required=TRUE, help="Output directory")
parser$add_argument('--nbases', type="numeric", default=1e+08, help="Number of bases to use for learning error model")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
# Set up logging
log_threshold(args$log_level)
log_appender(appender_console)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

if (!dir.exists(args$dout)) {
    dir.create(args$dout, recursive=TRUE)
}

# Learn errors for 1 or Merged using all loaded files
error_model <- learnErrors(args$filt_1, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand, nbases = args$nbases)

# If Reverse reads are not provided, save the error model as err_model.RDS
if (is.null(args$filt_2)) {
    saveRDS(error_model, file.path(args$dout, "err_model.RDS"))
} else {
    saveRDS(error_model, file.path(args$dout, "err_F_model.RDS"))
}

if (!is.null(args$filt_2)) {
    # Learn errors for 2 using all loaded files
    error_model <- learnErrors(args$filt_2, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand, nbases = args$nbases)
    saveRDS(error_model, file.path(args$dout, "err_R_model.RDS"))
}

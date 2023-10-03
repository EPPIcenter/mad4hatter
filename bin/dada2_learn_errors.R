library(logger)
log_threshold(WARN)
log_appender(appender_console)

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

load_library("lobstr")
# Define profiling function
profile_function <- function(func, ...) {
  # Record the memory size of the arguments
  mem_args <- lobstr::obj_size(...)

  # Record start time and memory
  time_start <- proc.time()[["elapsed"]]
  mem_start <- lobstr::mem_used()

  # Evaluate the function
  result <- func(...)

  # Record end time and memory
  time_end <- proc.time()[["elapsed"]]
  mem_end <- lobstr::mem_used()

  # Return results
  list(result = result,
      mem_args = mem_args,
      mem_diff = (mem_end - mem_start) - mem_args, 
      duration = time_end - time_start)
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
log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

if (!dir.exists(args$dout)) {
    dir.create(args$dout, recursive=TRUE)
}

# Learn errors for 1 or Merged using all loaded files
profiling_result_F <- profile_function(learnErrors, args$filt_1, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand, nbases = args$nbases)

log_info("Time taken for learning errors on forward reads: {profiling_result_F$duration}")
log_info("Memory change during learning errors on forward reads: {profiling_result_F$mem_diff}")

error_model <- profiling_result_F$result

# If Reverse reads are not provided, save the error model as err_model.RDS
output_file <- if (is.null(args$filt_2)) "err_model.RDS" else "err_F_model.RDS"
saveRDS(error_model, file.path(args$dout, output_file))

if (!is.null(args$filt_2)) {
    # Profile learnErrors for 2 using all loaded files
    profiling_result_R <- profile_function(learnErrors, args$filt_2, multithread=args$ncores, MAX_CONSIST=args$maxConsist, randomize=args$rand, nbases = args$nbases)
    
    log_info("Time taken for learning errors on reverse reads: {profiling_result_R$duration}")
    log_info("Memory change during learning errors on reverse reads: {profiling_result_R$mem_diff}")
    
    error_model <- profiling_result_R$result
    saveRDS(error_model, file.path(args$dout, "err_R_model.RDS"))
}

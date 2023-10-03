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
       mem_diff = (mem_end - mem_start) - mem_args, 
       duration = time_end - time_start)
}

load_library("argparse")
load_library("dada2")

# Argument Parsing and Setup
parser <- ArgumentParser(description='Perform Dereplication')
parser$add_argument('--filtFs-path', type="character", required=TRUE, help="Path to filtered F reads", nargs="+")
parser$add_argument('--filtRs-path', type="character", required=TRUE, help="Path to filtered R reads", nargs="+")
parser$add_argument('--filter-metadata', type="character", required=TRUE, help="Path to filter metadata")
parser$add_argument('--dout', type="character", required=TRUE, help="Output directory")
parser$add_argument('--ncores', type="numeric", required=TRUE, help="Number of threads to use")
parser$add_argument('--log-level', type="character", default = "INFO", help = "Log level. Default is INFO.")

args <- parser$parse_args()
# Set up logging
log_level_arg <- match.arg(args$log_level, c("DEBUG", "INFO", "WARN", "ERROR", "FATAL"))
log_threshold(log_level_arg)
args_string <- paste(sapply(names(args), function(name) {
  paste(name, ":", args[[name]])
}), collapse = ", ")

log_debug(paste("Arguments parsed successfully:", args_string))

# Load filter metadata
filter_metadata <- readRDS(args$filter_metadata)

if (!dir.exists(args$dout)) {
    dir.create(args$dout, recursive=TRUE)
}

# Profile the dereplication step
profiling_result <- profile_function(function(args) {
    # Dereplicate and save for each pair of filtFs and filtRs
    for(i in seq_along(args$filtFs_path)) {

        # Capture output and messages from derepFastq for F reads
        derepF_output <- capture.output({
            derepF <- derepFastq(args$filtFs_path[i], verbose = TRUE)
        }, type = "message")

        # If there's any message from derepFastq, log it
        if(length(derepF_output) > 0) {
          log_info(paste("Dereplication F reads [", args$filtFs_path[i], "]:", paste(derepF_output, collapse = "\n")))
        }
        # Construct unique filenames for saving F reads
        outFileF <- file.path(args$dout, paste0(basename(args$filtFs_path[i]), "_dereped.RDS"))
        saveRDS(derepF, file=outFileF)

        # Perform Dereplication for R reads
        derepR_output <- capture.output({
          derepR <- derepFastq(args$filtRs_path[i], verbose = TRUE)
        }, type = "message")

        if(length(derepR_output) > 0) {
          log_info(paste("Dereplication R reads [", args$filtRs_path[i], "]:", paste(derepR_output, collapse = "\n")))
        }
        # Construct unique filenames for saving R reads
        outFileR <- file.path(args$dout, paste0(basename(args$filtRs_path[i]), "_dereped.RDS"))
        saveRDS(derepR, file=outFileR)

    }
}, args)

# Log the profiling results
log_info("Total time taken for the loop: {profiling_result$duration}")
log_info("Memory change during the loop: {profiling_result$mem_diff}")
# Function To Push Commands Julia To Simulator

#' @title run_cmd
#' @description The run_cmd function takes one or more GraphML files and
#' simulation parameters, pushes them to the external Julia diffusion simulator,
#' and writes results to JSON and CSV files.
#' 
#' @param cmd File path to the executable.
#' @param mode A character value specifying the simulator mode. Options include:
#'   "single" for a single network, "sweep" for parameter sweeps, and "multi"
#'   for multiple networks.
#' @param graphml File path to a GraphML file or a comma-separated string of file
#'   paths for multiple networks.
#' @param infected Vector of seed node IDs.
#' @param inf_r Infection rate.
#' @param rec_t Recovery time (in time steps).
#' @param max_time Maximum simulation duration.
#' @param p_symp Probability of being symptomatic once infected.
#' @param b_int Baseline probability of activity.
#' @param b_close Edge weight effect on activity probability.
#' @param b_cxn_peers Effect of peers’ infection status on activity probability.
#' @param b_cxn_total Effect of global infection prevalence on activity probability.
#' @param b_cxn_symp Effect of being symptomatic on activity probability.
#' @param b_cls_x_smp Interaction of being symptomatic and edge weight on activity probability.
#' @param num_iter Number of simulation runs (for sweep mode).
#' @param transmission Transmission method. Use "weighted". 
#' @param seed Integer seed for reproducibility.
#' @param out File path for the JSON metadata output.
#' @param sas_transformation Logical; whether to perform SAS-style network transformation
#'   before simulation (default = FALSE).
#' @param echo Logical; whether to print the simulation command and output to the console
#'   (default = TRUE).
#' @param fail_on_status Logical; whether to stop execution if the simulator exits with a
#'   non-zero status (default = TRUE).
#'
#' @details
#' The function executes a call to an external Julia diffusion
#' simulator and returns two output files:
#' \itemize{
#'   \item a JSON metadata file (specified by \code{out}) summarizing the
#'     simulation results and referencing the infection log, and
#'   \item a CSV infection log file containing time-step results.
#' }
#'
#' The JSON file includes overall simulation metadata such as:
#' \itemize{
#'   \item \code{n_susceptible} – number of nodes remaining susceptible at simulation end.
#'   \item \code{n_infected} – number of nodes infected at simulation end.
#'   \item \code{n_recovered} – number of nodes recovered at simulation end.
#'   \item \code{infection_log} – file path to the infection log CSV.
#' }
#'
#' The infection log CSV records the state of the simulation at each time step:
#' \itemize{
#'   \item \code{time} – simulation time step.
#'   \item \code{n_infected} – number of nodes infected at each time step.
#'   \item \code{prop_ever_infected} – proportion of nodes ever infected.
#'   \item \code{prop_currently_infected} – proportion of nodes currently infected.
#'   \item \code{n_recovered} – number of nodes recovered.
#'   \item \code{prop_recovered} – proportion of nodes recovered.
#' }
#'
#' Users can access the infection log path directly from the JSON metadata
#' using \code{jsonlite::read_json}, or construct the file name manually:
#' \preformatted{
#' single_csv <- sub("\\\\.json$", "", SINGLE_OUT)
#' single_csv <- paste0(single_csv, "_infection_log.csv")
#' results <- readr::read_csv(single_csv)
#' }
#'
#' The two files together provide complete simulation output: the JSON
#' summarizes global outcomes, and the CSV contains detailed time-step data.
#'
#' @return
#' A character vector containing the simulator’s standard output,
#' with the following attributes:
#' \describe{
#'   \item{\code{status}}{Integer exit status from the simulator (0 = success).}
#'   \item{\code{stderr}}{Character vector containing any simulator error messages.}
#' }
#'
#' The simulator also writes two files to disk:
#' \itemize{
#'   \item a JSON metadata file (specified by \code{out}), which records
#'     simulation parameters and summary statistics, and
#'   \item a CSV infection log (referenced in the JSON metadata) that contains
#'     time-step infection and recovery data.
#' }
#'
#' The JSON and CSV files are not returned directly by the function but can be
#' read from disk after execution using \code{jsonlite::read_json} and
#' \code{readr::read_csv}, respectively.
#'
#' @examples
#' \dontrun{
#' ## Example 1: Single-network run
#' SINGLE_OUT <- tempfile(fileext = ".json")
#'
#' run_cmd(
#'   cmd = "diffusion_sim.exe",
#'   mode = "single",
#'   graphml = "example_network.graphml",
#'   infected = c(1, 2, 3),
#'   inf_r = 0.02,
#'   rec_t = 14,
#'   max_time = 200,
#'   p_symp = 0.75,
#'   b_int = -0.1,
#'   b_close = 1.0,
#'   b_cxn_peers = -0.5,
#'   b_cxn_total = -3.5,
#'   b_cxn_symp = -1.5,
#'   b_cls_x_smp = -0.1,
#'   transmission = "weighted",
#'   seed = 12345,
#'   out = SINGLE_OUT
#' )
#'
#' # Read infection log CSV and metadata JSON
#' single_csv <- sub("\\.json$", "", SINGLE_OUT)
#' single_csv <- paste0(single_csv, "_infection_log.csv")
#' results_single <- readr::read_csv(single_csv, show_col_types = FALSE)
#' meta_single <- jsonlite::read_json(SINGLE_OUT, simplifyVector = TRUE)
#'
#' ## Example 2: Parameter sweep run
#' SWEEP_OUT <- tempfile(fileext = ".json")
#'
#' run_cmd(
#'   cmd = "diffusion_sim.exe",
#'   mode = "sweep",
#'   graphml = "example_network.graphml",
#'   infected = 1,
#'   inf_r = 0.02,
#'   rec_t = 14,
#'   max_time = 200,
#'   p_symp = 0.75,
#'   b_int = -0.1,
#'   b_close = 1.0,
#'   b_cxn_peers = "-0.5,-0.75",
#'   b_cxn_total = "-3.5,-5.0",
#'   b_cxn_symp = "-1.5,-2.0",
#'   b_cls_x_smp = -0.1,
#'   num_iter = 2,
#'   transmission = "weighted",
#'   seed = 12345,
#'   out = SWEEP_OUT
#' )
#'
#' sweep_csv <- sub("\\.json$", "", SWEEP_OUT)
#' sweep_csv <- paste0(sweep_csv, "_infection_log.csv")
#' results_sweep <- readr::read_csv(sweep_csv, show_col_types = FALSE)
#' meta_sweep <- jsonlite::read_json(SWEEP_OUT, simplifyVector = TRUE)
#'
#' ## Example 3: Multi-network run
#' MULTI_OUT <- tempfile(fileext = ".json")
#' graphml_files <- paste(
#'   c("network1.graphml", "network2.graphml", "network3.graphml"),
#'   collapse = ","
#' )
#'
#' run_cmd(
#'   cmd = "diffusion_sim.exe",
#'   mode = "multi",
#'   graphml = graphml_files,
#'   infected = c(1, 3, 9),
#'   inf_r = 0.02,
#'   rec_t = 14,
#'   max_time = 200,
#'   p_symp = 0.75,
#'   b_int = -0.1,
#'   b_close = 1.0,
#'   b_cxn_peers = -0.5,
#'   b_cxn_total = -3.5,
#'   b_cxn_symp = -1.5,
#'   b_cls_x_smp = -0.1,
#'   transmission = "weighted",
#'   seed = 12345,
#'   out = MULTI_OUT
#' )
#'
#' multi_csv <- sub("\\.json$", "", MULTI_OUT)
#' multi_csv <- paste0(multi_csv, "_infection_log.csv")
#' results_multi <- readr::read_csv(multi_csv, show_col_types = FALSE)
#' meta_multi <- jsonlite::read_json(MULTI_OUT, simplifyVector = TRUE)
#' }
#'
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a
assert_file <- function(path, msg = NULL) {
  if (!file.exists(path)) stop(msg %||% paste0("Expected file not found: ", path))
}
if (!requireNamespace("readr", quietly = TRUE)) {
      if (isTRUE(auto_install)) {
        message("Installing 'readr'…")
        install.packages("readr")
      }
}
run_cmd <- function(cmd, mode, graphml, infected, inf_r, rec_t,
                    max_time, p_symp, b_int, b_close, b_cxn_peers,
                    b_cxn_total, b_cxn_symp, b_cls_x_smp, transmission,
                    seed, out, num_iter = NULL, sas_transformation = FALSE,
                    echo = TRUE, fail_on_status = TRUE) {
    #   Helper to collapse vectors for CLI
        collapse_vec <- function(x) paste(x, collapse = ",")
    
    #   Build CLI arguments
        args <- c(
          mode,
          "--graphml", shQuote(graphml),
          c("--infected", as.character(infected)),
          "--inf-r", as.character(inf_r),
          "--rec-t", as.character(rec_t),
          "--max-time", as.character(max_time),
          "--p-symp", as.character(p_symp),
          "--b-int", as.character(b_int),
          "--b-close", as.character(b_close),
          "--b-cxn-peers", collapse_vec(b_cxn_peers),
          "--b-cxn-total", collapse_vec(b_cxn_total),
          "--b-cxn-symp", collapse_vec(b_cxn_symp),
          "--b-cls-x-smp", as.character(b_cls_x_smp),
          "--transmission", transmission,
          "--seed", as.character(seed),
          "--out", out
        )
        
    #   Optional arguments
        if (!is.null(num_iter)) args <- c(args, "--num-iter", as.character(num_iter))
        if (isTRUE(sas_transformation)) args <- c(args, "--sas-transformation")
      
    #   Pretty echo of the command
        if (echo) cat("\n$ ", shQuote(cmd), " ", paste(shQuote(args), collapse = " "), "\n", sep = "")
    
    #   Resolve from PATH if needed
        resolved <- Sys.which(cmd)
        if (nzchar(resolved)) cmd <- resolved
    
    #   Capture both stdout and stderr to files, then read back
        outf <- tempfile("stdout_"); errf <- tempfile("stderr_")
        on.exit({ unlink(c(outf, errf), force = TRUE) }, add = TRUE)
        
        status <- suppressWarnings(
          system2(cmd, args, stdout = outf, stderr = errf, wait = TRUE)
        )
        out <- if (file.exists(outf)) readLines(outf, warn = FALSE) else character()
        err <- if (file.exists(errf)) readLines(errf, warn = FALSE) else character()
    
    #   Echo output for visibility
        if (echo && length(out)) cat(paste(out, collapse = "\n"), "\n", sep = "")
        if (echo && length(err)) cat(paste(err, collapse = "\n"), "\n", sep = "")
    
        if (is.na(status)) status <- 0L
        if (fail_on_status && status != 0L) {
          stop(
            "Non-zero exit status (", status, ") for: ", cmd, " ", paste(args, collapse = " "),
            if (length(err)) paste0("\nStderr:\n", paste(err, collapse = "\n")) else ""
          )
    }
  
  #   Return stdout, with useful attrs
  structure(out, status = status, stderr = err)
}
    
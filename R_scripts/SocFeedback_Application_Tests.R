#   Test Script for the SocFeedback Simulator
#   Jonathan H. Morgan, Ph.D.
#   3 October 2025

#   Clear Out Console Script
    cat("\014")
    rm(list = ls(all.names = TRUE))

#   Options
    options(stringsAsFactors = FALSE)
    options(mc.cores = parallel::detectCores())

#   Set Working Directory
    setwd("/workspace/data/sim_test_data")
    getwd()

#   Pointing to the Simulator Executable
    SIM_BIN <- "/workspace/ergm-based-diffusion-simulations/julia_env/build/diffusion_sim_app/bin/diffusion_sim"
    if (!file.exists(SIM_BIN)) SIM_BIN <- "diffusion_sim"  # fallback to PATH

#   Pointing to Data Directories for Example Cases
    SIM_BIN   <- "/workspace/ergm-based-diffusion-simulations/julia_env/build/diffusion_sim_app/bin/diffusion_sim"
    DATA_DIR   <- "/workspace/data/sim_test_data"
    GRAPHML    <- file.path(DATA_DIR, "WeakCore1_3.2_2.graphml")
    SINGLE_OUT <- file.path(DATA_DIR, "result_single.json")
    SWEEP_OUT  <- "/workspace/ergm-based-diffusion-simulations/result_sweep.json"
    MULTI_DIR  <- file.path(DATA_DIR, "synthetic_village_m_graphml")
    MULTI_OUT  <- "/workspace/ergm-based-diffusion-simulations/multi_test_results.json"
    SAS_CSV   <- "/workspace/data/sim_test_data/SAS_Simulation_Outputs.csv"
    STOP_IF_MISS <- FALSE

#################
#   FUNCTIONS   #
#################

#   Command Helper Function
    `%||%` <- function(a, b) if (is.null(a)) b else a

    assert_file <- function(path, msg = NULL) {
    if (!file.exists(path)) stop(msg %||% paste0("Expected file not found: ", path))
    }

#   Checking if readr is present
    if (!requireNamespace("readr", quietly = TRUE)) {
        if (isTRUE(auto_install)) {
            message("Installing 'readr'…")
            install.packages("readr")
        }
    }

#   Function for Pushing Commands to the Simulator
    run_cmd <- function(cmd, args = character(), echo = TRUE, fail_on_status = TRUE) {
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

#   Build the SAS Trial Design Matrix
    build_sas_design_vectors <- function(bproot = -0.5, bgroot = -3.5, bsroot = -1.5) {
        peer <- c(0.0, bproot, bproot - 0.25, bproot - 0.50)
        glob <- c(0.0, bgroot, bgroot - 1.5, bgroot - 3.0)
        self <- c(0.0, bsroot, bsroot - 0.5, bsroot - 1.0, bsroot - 1.5, bsroot - 2.0)
        list(peer = peer, glob = glob, self = self)
    }

#   SAS Replication Sweep
    run_sas_style_sweep <- function(sim_bin, graphml, sas_transformation = TRUE, infected = 1, inf_r = 0.02, rec_t = 14,
                                    max_time = 200, p_symp = 0.75, b_int = -0.1, b_close = 1.0,
                                    b_cls_x_smp = -0.1, design = build_sas_design_vectors(),
                                    transmission = "weighted", num_iter = 200, seed = 12345,
                                    out_json = tempfile(fileext = ".json"), echo_cmd = TRUE) {
        #   Ensure GraphML exists
            stopifnot(file.exists(graphml))
        
        #   Helper to collapse vectors for CLI ---
            v2s <- function(x) paste(x, collapse = ",")
        
        #   Build CLI arguments ---
          args <- c(
            "sweep",
            "--graphml", graphml,
            "--infected", as.character(infected),
            "--inf-r", as.character(inf_r),
            "--rec-t", as.character(rec_t),
            "--max-time", as.character(max_time),
            "--p-symp", as.character(p_symp),
            "--b-int", as.character(b_int),
            "--b-close", as.character(b_close),
            "--b-cxn-peers", v2s(design$peer),
            "--b-cxn-total", v2s(design$glob),
            "--b-cxn-symp", v2s(design$self),
            "--b-cls-x-smp", as.character(b_cls_x_smp),
            "--num-iter", as.character(num_iter),
            "--transmission", transmission,
            "--seed", as.character(seed),
            "--out", out_json
      )
          
      #   Add SAS transformation flag if TRUE
          if (sas_transformation) {
            args <- c(args, "--sas-transformation")
          }
      
      # --- Run the simulator ---
      run_cmd(sim_bin, args, echo = echo_cmd)
      
      # --- Read output CSV ---
      csv_path <- sub("\\.json$", "_infection_log.csv", out_json)
      assert_file(csv_path, "Sweep infection_log CSV not found.")
      
      df <- readr::read_csv(csv_path, show_col_types = FALSE)
      
      # --- Return stacked DataFrame ---
      df
    }
    

#   Simulator Results Aggregation
    julia_medians_from_df <- function(df) {
        #   Expected columns from CLI: time, prop_ever_infected, prop_currently_infected, prop_recovered
            req <- c("time","b_cxn_peers","b_cxn_total","b_cxn_symp",
                    "prop_ever_infected","prop_currently_infected","prop_recovered")
            missing <- setdiff(req, names(df))
            if (length(missing)) stop("Missing columns in Julia (CLI) CSV: ", paste(missing, collapse = ", "))

        #   Rename to match the Julia doc naming used in comparisons
            d <- df
            names(d)[names(d) == "b_cxn_peers"]   <- "b_cxn_peer"
            names(d)[names(d) == "b_cxn_total"]   <- "b_cxn_global"
            names(d)[names(d) == "b_cxn_symp"]    <- "b_self"
            names(d)[names(d) == "prop_ever_infected"]      <- "propeverinf"
            names(d)[names(d) == "prop_currently_infected"] <- "propcurrinf"
            names(d)[names(d) == "prop_recovered"]          <- "proprec"

        #   Medians per (time, coef triple)
            agg_key <- c("time","b_cxn_peer","b_cxn_global","b_self")
            med_time <- aggregate(
                d[c("propeverinf","propcurrinf","proprec")],
                by = d[agg_key],
                FUN = stats::median, na.rm = TRUE
            )
            names(med_time)[names(med_time) == "propeverinf"] <- "propeverinf_med_julia"
            names(med_time)[names(med_time) == "propcurrinf"] <- "propcurrinf_med_julia"
            names(med_time)[names(med_time) == "proprec"]     <- "proprec_med_julia"

        #   Finals per (coef triple × iter): CLI doesn't always carry 'itter'; tolerate absence
            iter_col <- intersect(c("iter","itter"), names(d))
            by_keys  <- c("b_cxn_peer","b_cxn_global","b_self", iter_col)
            if (length(iter_col) == 0) {
                #   Fabricate an iteration id = 1 for all
                    d$itter <- 1L
                    by_keys <- c("b_cxn_peer","b_cxn_global","b_self","itter")
            }

            final_per_run <- do.call(rbind, lapply(split(d, d[by_keys], drop = TRUE), function(g) {
                g <- g[order(g$time), , drop = FALSE]
                g[nrow(g), c("b_cxn_peer","b_cxn_global","b_self", iter_col, "time",
                            "propeverinf","propcurrinf","proprec")]
            }))
            names(final_per_run)[names(final_per_run) == "time"]          <- "time_final"
            names(final_per_run)[names(final_per_run) == "propeverinf"]   <- "propeverinf_final"
            names(final_per_run)[names(final_per_run) == "propcurrinf"]   <- "propcurrinf_final"
            names(final_per_run)[names(final_per_run) == "proprec"]       <- "proprec_final"

        #   Medians of finals across iterations
            med_final <- aggregate(
                final_per_run[c("time_final","propeverinf_final","propcurrinf_final","proprec_final")],
                by = final_per_run[c("b_cxn_peer","b_cxn_global","b_self")],
                FUN = stats::median, na.rm = TRUE
            )

        #   Return Medians
            return(list(med_time = med_time, final_per_run = final_per_run, med_final = med_final))
    }

#   Read + aggregate SAS CSV (align names) ----------
    read_sas_and_aggregate <- function(sas_csv) {
        #   Pulling-In SAS Results
            df <- readr::read_csv(sas_csv)

        #   Expected columns in SAS CSV (per your doc)
            keep <- c("ITTER","BETA","TIME","PROPEVERINF","PROPCURRINF","PROPREC",
                        "B_INT","B_CXN_PEER","B_CXN_GLOBAL","B_CLOSE","B_SELF","B_CLS_X_SMP","SEEDNODE")
            missing <- setdiff(keep, names(df))
            if (length(missing)) stop("SAS CSV missing columns: ", paste(missing, collapse = ", "))

            d <- df[keep]
            names(d) <- tolower(names(d))
            names(d)[names(d) == "b_cxn_peer"]   <- "b_cxn_peer"
            names(d)[names(d) == "b_cxn_global"] <- "b_cxn_global"
            names(d)[names(d) == "b_self"]       <- "b_self"

        #   Per-time medians (by coef triple)
            agg_key <- c("time","b_cxn_peer","b_cxn_global","b_self")
            med_time <- aggregate(
                d[c("propeverinf","propcurrinf","proprec","beta","b_int","b_close","b_cls_x_smp")],
                by = d[agg_key],
                FUN = function(x) if (is.numeric(x)) stats::median(x, na.rm = TRUE) else x[1]
            )
            names(med_time)[names(med_time) == "propeverinf"] <- "propeverinf_med_sas"
            names(med_time)[names(med_time) == "propcurrinf"] <- "propcurrinf_med_sas"
            names(med_time)[names(med_time) == "proprec"]     <- "proprec_med_sas"

        #   Return SAS Aggregated Data
            return(list(df = d, med_time = med_time))
    }

#   Compare SAS vs Julia (medians joined on time & coefs)
    compare_sas_vs_julia <- function(med_time_sas, med_time_julia,
                                    constants = list(beta = 0.02, p_symp = 0.75,
                                                    b_int = -0.1, b_close = 1.0, b_cls_x_smp = -0.1,
                                                    rec_t = 14, maxtime = 200)) {
        #   Creating Comparison Table
            key <- c("time","b_cxn_peer","b_cxn_global","b_self")
            comp <- dplyr::left_join(med_time_sas, med_time_julia, by=key)

            comp$propeverinf_delta <- abs(comp$propeverinf_med_sas - comp$propeverinf_med_julia)
            comp$propcurrinf_delta <- abs(comp$propcurrinf_med_sas - comp$propcurrinf_med_julia)
            comp$proprec_delta     <- abs(comp$proprec_med_sas     - comp$proprec_med_julia)

        #   Append constants
            for (nm in names(constants)) comp[[nm]] <- constants[[nm]]

        #   Reorder columns: constants → coefs → time → medians → deltas
            comp <- comp[, c("beta","p_symp","b_int","b_close","b_cls_x_smp","rec_t","maxtime",
                            "b_cxn_peer","b_cxn_global","b_self","time",
                            "propeverinf_med_sas","propcurrinf_med_sas","proprec_med_sas",
                            "propeverinf_med_julia","propcurrinf_med_julia","proprec_med_julia",
                            "propeverinf_delta","propcurrinf_delta","proprec_delta")]

        #   Return Comparison Data
            comp
    }

#   Plotting Curves
    plot_sas_vs_julia_panels <- function(comp, outcome = c("propeverinf", "propcurrinf", "proprec")) {
    
        #   Specifying Outcome to Plot
            outcome <- match.arg(outcome)

        #   Column names to use
            col_j <- switch(outcome,
                            propeverinf = "propeverinf_med_julia",
                            propcurrinf = "propcurrinf_med_julia",
                            proprec     = "proprec_med_julia")
            col_s <- switch(outcome,
                            propeverinf = "propeverinf_med_sas",
                            propcurrinf = "propcurrinf_med_sas",
                            proprec     = "proprec_med_sas")
            y_label <- switch(outcome,
                                propeverinf = "Proportion Ever Infected",
                                propcurrinf = "Proportion Currently Infected",
                                proprec     = "Proportion Recovered")

        #   Coerce types (safe)
            comp <- within(comp, {
                b_cxn_peer   <- as.numeric(b_cxn_peer)
                b_cxn_global <- as.numeric(b_cxn_global)
                b_self       <- as.numeric(b_self)
                time         <- as.numeric(time)
            })

            peer_levels <- sort(unique(comp$b_cxn_peer))
            glob_levels <- sort(unique(comp$b_cxn_global))
            self_levels <- sort(unique(comp$b_self))

        #   Ticks helper (matches your Hmisc usage; installs if needed)
            dataplot_tick_function <- function(major_tick_length = 0.035, minor_tick_ratio = 0.25) {
                if (!requireNamespace("Hmisc", quietly = TRUE)) {
                install.packages("Hmisc", repos = "https://cloud.r-project.org")
                }
                Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = minor_tick_ratio)
                Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = -minor_tick_ratio)
                axis(2, tck = 1, tck = -major_tick_length, labels = FALSE)
                axis(1, tck = 1, tck = -major_tick_length, labels = FALSE)
            }

        #   Palette exactly like your base set, extended if needed
            pal_base <- c("#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#4477AA", "#BBBBBB")
            pal <- if (length(self_levels) <= length(pal_base)) {
                pal_base[seq_along(self_levels)]
            } else {
                grDevices::colorRampPalette(pal_base)(length(self_levels))
            }

        #   Laying-Out the Plot
            layout(matrix(1:(length(glob_levels) * length(peer_levels)),
                            nrow = length(glob_levels), byrow = TRUE))
            par(family = "serif", mar = c(3.5, 4, 2.5, 1.5), las = 1)

        #   Plotting
            for (gi in seq_along(glob_levels)) {
                #   Creating Plotting Variables
                    gval <- glob_levels[gi]

                #   Plotting Panels
                    for (pi in seq_along(peer_levels)) {
                        pval <- peer_levels[pi]
                        sub  <- comp[comp$b_cxn_global == gval & comp$b_cxn_peer == pval, , drop = FALSE]

                        #   y-range across SAS & Julia series
                            ylim <- range(c(sub[[col_j]], sub[[col_s]]), na.rm = TRUE)
                            ylim[1] <- min(0, ylim[1])

                        #   Creating Base Plot
                            plot(NA, type = "n",
                                xlim = range(sub$time, na.rm = TRUE), ylim = ylim,
                                xlab = " ", ylab = y_label, tck = 0.015, xaxt = 'n', bty = 'L', las = 1,
                                main = paste0("Peer effect = ", sprintf("%.2f", pval),
                                                "   |   Global effect = ", sprintf("%.2f", gval)))
                            mtext(side = 1, text = "Time", col = "black", line = 2.45, cex = 0.75, family = 'serif')
                            axis(1, padj = 0.75, tck = 0.015)
                            dataplot_tick_function()

                        #   One line per symptomatic level; Julia solid, SAS dashed
                            for (si in seq_along(self_levels)) {
                                sval <- self_levels[si]
                                ss <- sub[sub$b_self == sval, , drop = FALSE]
                                ss <- ss[order(ss$time), , drop = FALSE]
                                if (nrow(ss) == 0) next
                                lines(ss$time, ss[[col_j]], lwd = 2, lty = 1, col = pal[si])
                                lines(ss$time, ss[[col_s]], lwd = 2, lty = 2, col = pal[si])
                            }
                    }
            }

    }

#############
#   TESTS   #
#############
    
#   Graphml Test
    pajek_file <- file.path(DATA_DIR, "WeakCore1_3.2_2.net")
    graphml_file <- file.path(DATA_DIR, "test.graphml")
    g <- igraph::read_graph(pajek_file, format = "pajek")
    
    igraph::write_graph(g, graphml_file, format = "graphml")
    g_2 <- igraph::read_graph(pajek_file, format = "pajek")
    summary(g_2)

#   Test 0: Location Test
    cat("==> Using SIM_BIN: ", SIM_BIN, "\n", sep = "")
    if (identical(SIM_BIN, "diffusion_sim")) {
        #   confirm it's on PATH
            which_out <- suppressWarnings(system2("which", "diffusion_sim", stdout = TRUE))
            if (length(which_out) == 0L) {
                stop("diffusion_sim not found on PATH. Set SIM_BIN to the full path, or export PATH accordingly.")
            }
            cat("    Resolved on PATH: ", which_out, "\n", sep = "")
    }

#   Test 1: Help Test
    cat("\n== Test 1: CLI help ==\n")
    run_cmd(SIM_BIN, c("single", "--help"))
    run_cmd(SIM_BIN, c("sweep",  "--help"))
    run_cmd(SIM_BIN, c("multi",  "--help"))
    cat("✓ Help screens printed.\n")

#   Test 2: Single Network Test
    cat("\n== Test 2: single run ==\n")
    dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)
    stopifnot(file.exists(GRAPHML))

    # "--sas-transformation",     # uncomment if you want SAS Network Transformation
    single_args <- c("single", "--graphml", GRAPHML, "--infected", "1", "3", "9", "--inf-r", "0.02",
                     "--rec-t", "14", "--max-time", "200", "--p-symp", "0.75", "--b-int", "-0.1", "--b-close", "1.0",
                     "--b-cxn-peers", "-0.5", "--b-cxn-total", "-3.5", "--b-cxn-symp", "-1.5", "--b-cls-x-smp", "-0.1",
                     "--transmission", "weighted", "--seed", "12345", "--out", SINGLE_OUT)
    run_cmd(SIM_BIN, single_args)

    single_csv <- sub("\\.json$", "", SINGLE_OUT)
    single_csv <- paste0(single_csv, "_infection_log.csv")
    assert_file(single_csv, "Single-run CSV not found.")

    df_single <- readr::read_csv(single_csv)
    cat("Single rows:", nrow(df_single), "  cols:", paste(colnames(df_single), collapse = ", "), "\n")
    print(utils::head(df_single, 10))

#   TEST 3: Basic Sweep Run (comma-separated vectors) 
    cat("\n== Test 3: sweep run ==\n")
    sweep_args <- c(
    "sweep",
    "--graphml", GRAPHML,
    # "--sas-transformation",
    "--infected", "1",
    "--inf-r", "0.02",
    "--rec-t", "14",
    "--max-time", "200",
    "--p-symp", "0.75",
    "--b-int", "-0.1",
    "--b-close", "1.0",
    "--b-cxn-peers", "-0.5,-0.75",
    "--b-cxn-total", "-3.5,-5.0",
    "--b-cxn-symp", "-1.5,-2.0",
    "--b-cls-x-smp", "-0.1",
    "--num-iter", "2",
    "--transmission", "weighted",
    "--seed", "12345",
    "--out", SWEEP_OUT
    )
    run_cmd(SIM_BIN, sweep_args)

    sweep_csv <- sub("\\.json$", "", SWEEP_OUT)
    sweep_csv <- paste0(sweep_csv, "_infection_log.csv")
    assert_file(sweep_csv, "Sweep CSV not found.")
    df_sweep <- readr::read_csv(sweep_csv)
    head(df_sweep)
    nrow(df_sweep)

#   Test 4: Basic Aggregation Test
    cat("\n== Test 4: columns & grouping sanity ==\n")
    expected <- c("time","n_infected","prop_ever_infected","prop_currently_infected",
                "n_recovered","prop_recovered","b_cxn_peers","b_cxn_total","b_cxn_symp","iter")
    present  <- colnames(df_sweep)
    cat("Expected:\n  ", paste(expected, collapse = ", "), "\n", sep = "")
    cat("Present:\n  ", paste(present,  collapse = ", "), "\n", sep = "")
    missing_cols <- setdiff(expected, present)
    if (length(missing_cols)) stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))

#   Each (bp,bg,bs,iter) should produce a full trajectory
    traj_lengths <- aggregate(time ~ b_cxn_peers + b_cxn_total + b_cxn_symp + iter,
                              data = df_sweep, FUN = length)
    cat("distinct trajectories: ", nrow(traj_lengths), "\n", sep = "")
    print(utils::head(traj_lengths, 8))

#   Check each group’s min time == 0
    mins <- aggregate(time ~ b_cxn_peers + b_cxn_total + b_cxn_symp + iter,
                    data = df_sweep, FUN = min)
    if (!all(mins$time == 0)) stop("Some trajectories do not start at time==0")

#   Check that no time exceeds your horizon (e.g., 200)
    if (max(df_sweep$time, na.rm = TRUE) > 200) stop("Found time > maxtime")
    cat("✓ Columns present, grouping sanity checks passed.\n")

#   TEST 5: Processing Multiple Networks
    cat("\n== Test 5: multi over 10 GraphML files ==\n")
    stopifnot(dir.exists(MULTI_DIR))
    all_files <- list.files(MULTI_DIR, full.names = TRUE)
    graphml_files <- all_files[grepl("\\.graphml$", all_files, ignore.case = TRUE)]
    if (length(graphml_files) < 10) {
    stop("Expected at least 10 .graphml files in: ", MULTI_DIR)
    }
    graphml_files <- graphml_files[1:10]
    file_list <- paste(graphml_files, collapse = ",")

    multi_args <- c(
    "multi",
    "--graphml", file_list,
    "--infected", "1", "3", "9",
    "--inf-r", "0.02",
    "--rec-t", "14",
    "--max-time", "200",
    "--p-symp", "0.75",
    "--b-int", "-0.1",
    "--b-close", "1.0",
    "--b-cxn-peers", "-0.5",
    "--b-cxn-total", "-3.5",
    "--b-cxn-symp", "-1.5",
    "--b-cls-x-smp", "-0.1",
    "--transmission", "weighted",
    "--seed", "12345",
    "--out", MULTI_OUT
    )
    run_cmd(SIM_BIN, multi_args)
    multi_csv <- paste0(sub("\\.json$", "", MULTI_OUT), "_infection_log.csv")
    assert_file(multi_csv, "Multi-run CSV not found.")
    cat("✓ Multi-mode produced: ", multi_csv, "\n", sep = "")

    df_multi <- readr::read_csv(multi_csv)
    head(df_multi)
    nrow(df_multi)

#######################
#   SAS REPLICATION   #
#######################

#   Sanity Checks
    if (!file.exists(SIM_BIN)) {
        msg <- paste0("SIM_BIN not found: ", SIM_BIN)
    if (STOP_IF_MISS) stop(msg) else warning(msg)
    }
    if (!file.exists(GRAPHML)) {
        msg <- paste0("GRAPHML not found: ", GRAPHML)
    if (STOP_IF_MISS) stop(msg) else warning(msg)
    }

#   Build the SAS design (defaults: -0.5, -3.5, -1.5)
    design <- build_sas_design_vectors()

#   Run SAS Replication
    df_cli <- run_sas_style_sweep(
      sim_bin   = SIM_BIN,
      graphml   = GRAPHML,
      sas_transformation = TRUE,
      infected  = 1,                # single seed node per iteration
      inf_r     = 0.02,
      rec_t     = 14,
      max_time  = 200,
      p_symp    = 0.75,
      b_int     = -0.1,
      b_close   = 1.0,
      b_cls_x_smp = -0.1,
      design    = design,
      transmission = "weighted",
      seed      = 12345,
      out_json  = tempfile(fileext = ".json"),
      echo_cmd  = TRUE
    )
    
#   Aggregate Julia medians/finals
    jl <- julia_medians_from_df(df_cli)
    med_time_julia <- jl$med_time

#   Read + aggregate SAS CSV
    sas <- read_sas_and_aggregate(SAS_CSV)
    med_time_sas <- sas$med_time

#   Build comparison (adds constants for traceability)
    comp <- compare_sas_vs_julia(
    med_time_sas,
    med_time_julia,
    constants = list(beta = 0.02, p_symp = 0.75,
                    b_int = -0.1, b_close = 1.0, b_cls_x_smp = -0.1,
                    rec_t = 14, maxtime = 200)
    )

#   Quick look & optional plot
    cat("Comparison rows: ", nrow(comp), "\n")
    print(utils::head(comp, 12))

#   Plotting
    plot_sas_vs_julia_panels(comp, outcome="propeverinf")
#' @title parameter_sweep_plot
#' @description
#' The function generates faceted panels of infection curves across
#' combinations of peer, global, and symptomatic infection parameters.
#' Each solid line represents a distinct symptomatic parameter level.
#'
#' @param df Data frame containing sweep simulation results, typically
#' produced by the simulator in sweep mode.
#' Must include columns \code{time}, \code{b_cxn_peers}, \code{b_cxn_total},
#' \code{b_cxn_symp}, and one or more outcome variables.
#' @param outcome Character; one of \code{"prop_ever_infected"},
#' \code{"prop_currently_infected"}, or \code{"prop_recovered"}.
#'
#' @details
#' If the input includes multiple iterations (\code{iter} column present),
#' the function automatically computes the median value of each outcome
#' across replicates for each parameter combination and time step.
#'
#' @examples
#' parameter_sweep_plot(df_sweep, outcome = "prop_ever_infected")
#'
#' @export
parameter_sweep_plot <- function(df,
                                 outcome = c("prop_ever_infected",
                                             "prop_currently_infected",
                                             "prop_recovered")) {
  #	"""
  #	Args:
  #		df: data frame with sweep results including time, b_cxn_peers, b_cxn_total, b_cxn_symp.
  #		outcome: which infection proportion to plot ("prop_ever_infected", "prop_currently_infected", or "prop_recovered").
  #	Returns:
  #		No return. Produces faceted base R plots of Julia sweep results.
  #	Notes:
  #		Automatically aggregates by median if multiple iterations (iter column) detected.
  #	"""
  
  #	Validation
    outcome <- match.arg(outcome)
    required_cols <- c("time", "b_cxn_peers", "b_cxn_total", "b_cxn_symp")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols)) {
      stop("Missing expected columns: ", paste(missing_cols, collapse = ", "))
    }
  
  #	Aggregation
    if ("iter" %in% names(df)) {
      message("Detected multiple iterations; summarizing by median across replicates...")
      df <- df |>
        dplyr::group_by(b_cxn_peers, b_cxn_total, b_cxn_symp, time) |>
        dplyr::summarise(
          prop_ever_infected      = median(prop_ever_infected, na.rm = TRUE),
          prop_currently_infected = median(prop_currently_infected, na.rm = TRUE),
          prop_recovered          = median(prop_recovered, na.rm = TRUE),
          .groups = "drop"
        )
    }
  
  #	Y-Axis Label
    y_label <- switch(outcome,
                      prop_ever_infected      = "Proportion Ever Infected",
                      prop_currently_infected = "Proportion Currently Infected",
                      prop_recovered          = "Proportion Recovered")
  
  #	Type Coercion
    df <- within(df, {
      b_cxn_peers <- as.numeric(b_cxn_peers)
      b_cxn_total <- as.numeric(b_cxn_total)
      b_cxn_symp  <- as.numeric(b_cxn_symp)
      time        <- as.numeric(time)
    })
  
  #	Parameter Levels
    peer_levels <- sort(unique(df$b_cxn_peers))
    glob_levels <- sort(unique(df$b_cxn_total))
    self_levels <- sort(unique(df$b_cxn_symp))
  
  #	Tick Helper
    dataplot_tick_function <- function(major_tick_length = 0.035,
                                       minor_tick_ratio  = 0.25) {
      if (!requireNamespace("Hmisc", quietly = TRUE)) {
        install.packages("Hmisc", repos = "https://cloud.r-project.org")
      }
      Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio =  minor_tick_ratio)
      Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = -minor_tick_ratio)
      axis(2, tck = 1, tck = -major_tick_length, labels = FALSE)
      axis(1, tck = 1, tck = -major_tick_length, labels = FALSE)
    }
  
  #	Color Palette
    pal_base <- c("#66CCEE", "#228833", "#CCBB44", "#EE6677",
                  "#AA3377", "#4477AA", "#BBBBBB")
    pal <- if (length(self_levels) <= length(pal_base)) {
      pal_base[seq_along(self_levels)]
    } else {
      grDevices::colorRampPalette(pal_base)(length(self_levels))
    }
  
  #	Layout
    layout(matrix(1:(length(glob_levels) * length(peer_levels)),
                  nrow = length(glob_levels), byrow = TRUE))
    par(family = "serif", mar = c(3.5, 4, 2.5, 1.5), las = 1)
  
  #	Plot Panels
    for (gi in seq_along(glob_levels)) {
      gval <- glob_levels[gi]
      
      for (pi in seq_along(peer_levels)) {
        pval <- peer_levels[pi]
        sub <- df[df$b_cxn_total == gval &
                    df$b_cxn_peers == pval, , drop = FALSE]
        
        ylim <- range(sub[[outcome]], na.rm = TRUE)
        ylim[1] <- min(0, ylim[1])
        
        plot(NA, type = "n",
             xlim = range(sub$time, na.rm = TRUE), ylim = ylim,
             xlab = " ", ylab = y_label,
             tck = 0.015, xaxt = "n", bty = "L", las = 1,
             main = paste0("Peer effect = ", sprintf("%.2f", pval),
                           "   |   Global effect = ", sprintf("%.2f", gval)))
        
        mtext(side = 1, text = "Time", col = "black",
              line = 2.45, cex = 0.75, family = "serif")
        axis(1, padj = 0.75, tck = 0.015)
        dataplot_tick_function()
        
        #	Draw lines by symptomatic level
        for (si in seq_along(self_levels)) {
          sval <- self_levels[si]
          ss <- sub[sub$b_cxn_symp == sval, , drop = FALSE]
          ss <- ss[order(ss$time), , drop = FALSE]
          if (nrow(ss) == 0) next
          lines(ss$time, ss[[outcome]], lwd = 2, lty = 1, col = pal[si])
        }
      }
    }
  }

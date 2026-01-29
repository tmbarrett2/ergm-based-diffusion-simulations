#' @title plot_sir_single
#' @description
#' Produces a single-panel SIR (Susceptible-Infected-Recovered) infection curve
#' plot with loess-smoothed trajectories. If multiple iterations are present,
#' plots individual iteration curves in gray with a colored median trend line
#' for each compartment.
#'
#' @param df Data frame containing simulation results with columns:
#' \itemize{
#'   \item \code{time} - Time step (numeric)
#'   \item \code{prop_ever_infected} - Cumulative proportion infected (numeric)
#'   \item \code{prop_currently_infected} - Current proportion infected (numeric)
#'   \item \code{prop_recovered} - Proportion recovered (numeric)
#'   \item \code{iteration} - Iteration ID (optional, for multiple runs)
#' }
#'
#' @details
#' Susceptible proportion is calculated as \code{1 - prop_ever_infected}.
#' 
#' For data with an \code{iteration} column:
#' \itemize{
#'   \item Individual iteration trajectories shown as semi-transparent gray lines
#'   \item Median trajectory shown as bold colored line
#'   \item All curves smoothed using loess (span = 0.75 for iterations, 0.5 for median)
#' }
#' 
#' Colors: Blue = Susceptible, Red = Infected, Green = Recovered
#' 
#' Plot extends from time 0 to maximum time in data.
#' 
#' Requires \code{ggplot2}, \code{dplyr}, and \code{tidyr} packages.
#'
#' @return A ggplot object
#'
#' @examples
#' plot_sir_single(simulation_results)
#'
#' @export
	plot_sir_single <- function(df) {
		#	"""
		#	Args:
		#		df: data.frame with columns time, prop_ever_infected,
		#		    prop_currently_infected, prop_recovered, and optionally iteration
		#	Returns:
		#		ggplot2 object showing SIR trajectory
		#	Notes:
		#		If iteration column present, shows individual runs in gray with
		#		median trajectory in color. Otherwise shows single smoothed run.
		#	"""
		
		#	Validation
			required_cols <- c("time", "prop_ever_infected", 
			                   "prop_currently_infected", "prop_recovered")
			missing_cols <- setdiff(required_cols, names(df))
			if (length(missing_cols) > 0) {
				stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
			}
		
		#	Check for Required Packages
			if (!requireNamespace("ggplot2", quietly = TRUE)) {
				stop("Package 'ggplot2' is required")
			}
		
		#	Calculate Susceptible Proportion
			df <- df |>
				dplyr::mutate(
					time = as.numeric(time),
					prop_ever_infected = as.numeric(prop_ever_infected),
					prop_currently_infected = as.numeric(prop_currently_infected),
					prop_recovered = as.numeric(prop_recovered),
					prop_susceptible = 1 - prop_ever_infected
				)
		
		#	Define Time Sequence
			max_time <- max(df$time, na.rm = TRUE)
			time_seq <- seq(0, max_time, length.out = 200)
		
		#	Check for Multiple Iterations
			if ("iteration" %in% names(df)) {
				
				#	Fit Loess for Each Iteration
					iterations <- unique(df$iteration)
					pred_list <- lapply(iterations, function(iter) {
						df_iter <- df[df$iteration == iter, ]
						
						loess_s <- loess(prop_susceptible ~ time, data = df_iter, span = 0.75)
						loess_i <- loess(prop_currently_infected ~ time, data = df_iter, span = 0.75)
						loess_r <- loess(prop_recovered ~ time, data = df_iter, span = 0.75)
						
						data.frame(
							time = rep(time_seq, 3),
							iteration = iter,
							compartment = factor(
								rep(c("Susceptible", "Infected", "Recovered"), each = length(time_seq)),
								levels = c("Susceptible", "Infected", "Recovered")
							),
							proportion = c(
								predict(loess_s, newdata = data.frame(time = time_seq)),
								predict(loess_i, newdata = data.frame(time = time_seq)),
								predict(loess_r, newdata = data.frame(time = time_seq))
							)
						)
					})
				
				#	Clean Predictions
					df_pred <- dplyr::bind_rows(pred_list) |>
						dplyr::filter(!is.na(proportion)) |>
						dplyr::mutate(proportion = pmax(0, pmin(1, proportion)))
				
				#	Calculate Raw Median Trajectory
					df_median_raw <- df_pred |>
						dplyr::group_by(time, compartment) |>
						dplyr::summarise(proportion = median(proportion, na.rm = TRUE), .groups = "drop")
				
				#	Smooth Median Helper
					smooth_median <- function(comp) {
						med <- df_median_raw |> dplyr::filter(compartment == comp)
						loess_med <- loess(proportion ~ time, data = med, span = 0.5)
						pred <- predict(loess_med, newdata = data.frame(time = time_seq))
						
						#	Force to actual value at time zero
							t0_val <- med$proportion[which.min(abs(med$time - 0))]
							pred[1] <- t0_val
						
						data.frame(
							time = time_seq,
							compartment = factor(comp, levels = c("Susceptible", "Infected", "Recovered")),
							proportion = pred
						)
					}
				
				#	Build Smoothed Median
					df_median <- dplyr::bind_rows(
						smooth_median("Susceptible"),
						smooth_median("Infected"),
						smooth_median("Recovered")
					) |>
						dplyr::filter(!is.na(proportion)) |>
						dplyr::mutate(proportion = pmax(0, pmin(1, proportion)))
				
			} else {
				
				#	Single Run: Smooth the Data
					loess_s <- loess(prop_susceptible ~ time, data = df, span = 0.75)
					loess_i <- loess(prop_currently_infected ~ time, data = df, span = 0.75)
					loess_r <- loess(prop_recovered ~ time, data = df, span = 0.75)
					
					df_pred <- data.frame(
						time = rep(time_seq, 3),
						iteration = 1,
						compartment = factor(
							rep(c("Susceptible", "Infected", "Recovered"), each = length(time_seq)),
							levels = c("Susceptible", "Infected", "Recovered")
						),
						proportion = c(
							predict(loess_s, newdata = data.frame(time = time_seq)),
							predict(loess_i, newdata = data.frame(time = time_seq)),
							predict(loess_r, newdata = data.frame(time = time_seq))
						)
					)
					
					df_median <- df_pred |>
						dplyr::select(time, compartment, proportion)
			}
		
		#	Define Color Palette
			colors <- c("Susceptible" = "#4477AA", 
			            "Infected" = "#EE6677", 
			            "Recovered" = "#228833")
		
		#	Create Plot
			p <- ggplot2::ggplot(df_pred, ggplot2::aes(x = time, y = proportion, 
			                                           group = interaction(iteration, compartment))) +
				ggplot2::geom_line(linewidth = 0.5, alpha = 0.2, color = "gray50") +
				ggplot2::geom_line(
					data = df_median, 
					ggplot2::aes(x = time, y = proportion, color = compartment, group = compartment),
					linewidth = 1.2, 
					inherit.aes = FALSE
				) +
				ggplot2::scale_color_manual(values = colors) +
				ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
				ggplot2::scale_x_continuous(limits = c(0, max_time), expand = c(0, 0), oob = scales::oob_keep) +
				ggplot2::coord_cartesian(xlim = c(0, max_time * 1.02), clip = "off") +
				ggplot2::labs(x = "Time", y = "Proportion", color = NULL) +
				ggplot2::theme_classic(base_family = "serif", base_size = 24) +
				ggplot2::theme(
					legend.position = "top",
					legend.direction = "horizontal",
					axis.text = ggplot2::element_text(size = 18),
					axis.title = ggplot2::element_text(size = 24),
					legend.text = ggplot2::element_text(size = 24),
					axis.line = ggplot2::element_blank(),
					panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5)
				)
		
		return(p)
	}

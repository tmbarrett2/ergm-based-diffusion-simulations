#' @title plot_sir_multi
#' @description
#' Produces a multi-panel SIR plot with one panel per network using base R graphics.
#' Each panel shows loess-smoothed trajectories across multiple iterations.
#'
#' @param df Data frame containing simulation results with columns:
#' \itemize{
#'   \item \code{network_id} - Network identifier (character/factor)
#'   \item \code{time} - Time step (numeric)
#'   \item \code{prop_ever_infected} - Cumulative proportion infected (numeric)
#'   \item \code{prop_currently_infected} - Current proportion infected (numeric)
#'   \item \code{prop_recovered} - Proportion recovered (numeric)
#'   \item \code{iteration} - Iteration ID (numeric)
#' }
#'
#' @return NULL (produces base R plot)
#'
#' @examples
#' plot_sir_multi(simulation_results)
#'
#' @export
#	Multi-Network SIR Plot
	plot_sir_multi <- function(df) {
		#	"""
		#	Args:
		#		df: data.frame with columns network_id, time, prop_ever_infected,
		#		    prop_currently_infected, prop_recovered, iteration
		#	Returns:
		#		NULL (invisible); generates multi-panel plot as side effect
		#	Notes:
		#		Creates 4-column layout showing SIR trajectories for each network.
		#		Individual iterations shown in gray; median trajectory in color.
		#	"""
		
		#	Validation
			required_cols <- c("network_id", "time", "prop_ever_infected", 
			                   "prop_currently_infected", "prop_recovered", "iteration")
			missing_cols <- setdiff(required_cols, names(df))
			if (length(missing_cols) > 0) {
				stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
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
		
		#	Initialize Storage
			networks <- sort(unique(df$network_id))
			n_networks <- length(networks)
			all_pred <- list()
			all_median <- list()
		
		#	Process Each Network
			for (net in networks) {
				df_net <- df |> dplyr::filter(network_id == net)
				iterations <- unique(df_net$iteration)
				
				#	Fit Loess for Each Iteration
					pred_list <- lapply(iterations, function(iter) {
						df_iter <- df_net[df_net$iteration == iter, ]
						
						loess_s <- loess(prop_susceptible ~ time, data = df_iter, span = 0.75)
						loess_i <- loess(prop_currently_infected ~ time, data = df_iter, span = 0.75)
						loess_r <- loess(prop_recovered ~ time, data = df_iter, span = 0.75)
						
						data.frame(
							network_id = net,
							time = rep(time_seq, 3),
							iteration = iter,
							compartment = c(rep("Susceptible", length(time_seq)),
							                rep("Infected", length(time_seq)),
							                rep("Recovered", length(time_seq))),
							proportion = c(
								predict(loess_s, newdata = data.frame(time = time_seq)),
								predict(loess_i, newdata = data.frame(time = time_seq)),
								predict(loess_r, newdata = data.frame(time = time_seq))
							)
						)
					})
				
				#	Clean Predictions
					df_pred_net <- dplyr::bind_rows(pred_list) |>
						dplyr::filter(!is.na(proportion)) |>
						dplyr::mutate(proportion = pmax(0, pmin(1, proportion)))
					
					all_pred[[net]] <- df_pred_net
				
				#	Calculate Raw Median Trajectory
					df_median_raw <- df_pred_net |>
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
							network_id = net,
							time = time_seq,
							compartment = comp,
							proportion = pred
						)
					}
				
				#	Build Smoothed Median
					df_median_net <- dplyr::bind_rows(
						smooth_median("Susceptible"),
						smooth_median("Infected"),
						smooth_median("Recovered")
					) |>
						dplyr::filter(!is.na(proportion)) |>
						dplyr::mutate(proportion = pmax(0, pmin(1, proportion)))
					
					all_median[[net]] <- df_median_net
			}
		
		#	Combine All Networks
			df_pred_all <- dplyr::bind_rows(all_pred)
			df_median_all <- dplyr::bind_rows(all_median)
		
		#	Tick Helper Function
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
		
		#	Define Color Palette
			colors <- list(
				Susceptible = "#4477AA",
				Infected = "#EE6677",
				Recovered = "#228833"
			)
		
		#	Configure Layout
			n_cols <- 4
			n_rows <- ceiling(n_networks / n_cols)
			layout(matrix(1:(n_rows * n_cols), nrow = n_rows, byrow = TRUE))
			par(family = "serif", mar = c(3.5, 4, 2.5, 1.5), las = 1)
		
		#	Plot Each Network
			for (i in 1:n_networks) {
				net <- networks[i]
				
				#	Subset Data
					sub_pred <- df_pred_all[df_pred_all$network_id == net, , drop = FALSE]
					sub_median <- df_median_all[df_median_all$network_id == net, , drop = FALSE]
				
				#	Create Empty Plot
					plot(NA, type = "n",
					     xlim = c(0, max_time), ylim = c(0, 1),
					     xlab = " ", ylab = "Proportion",
					     tck = 0.015, xaxt = "n", bty = "L", las = 1,
					     main = net)
					
					mtext(side = 1, text = "Time", col = "black",
					      line = 2.45, cex = 0.75, family = "serif")
					axis(1, padj = 0.75, tck = 0.015)
					dataplot_tick_function()
				
				#	Plot Individual Trajectories
					for (iter in unique(sub_pred$iteration)) {
						for (comp in c("Susceptible", "Infected", "Recovered")) {
							ss <- sub_pred[sub_pred$iteration == iter & sub_pred$compartment == comp, , drop = FALSE]
							ss <- ss[order(ss$time), , drop = FALSE]
							if (nrow(ss) == 0) next
							lines(ss$time, ss$proportion, lwd = 0.5, lty = 1, col = rgb(0.5, 0.5, 0.5, 0.2))
						}
					}
				
				#	Plot Median Trajectories
					for (comp in c("Susceptible", "Infected", "Recovered")) {
						ss <- sub_median[sub_median$compartment == comp, , drop = FALSE]
						ss <- ss[order(ss$time), , drop = FALSE]
						if (nrow(ss) == 0) next
						lines(ss$time, ss$proportion, lwd = 2, lty = 1, col = colors[[comp]])
					}
			}
		
		#	Fill Remaining Panels
			n_empty <- (n_rows * n_cols) - n_networks
			if (n_empty > 0) {
				for (i in 1:n_empty) {
					plot.new()
				}
			}
		
		#	Reset Layout
			layout(1)
			par(mar = c(5.1, 4.1, 4.1, 2.1))
		
		invisible(NULL)
	}

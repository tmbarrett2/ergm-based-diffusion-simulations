#SIR Social Feedback Simulation
#Jonathan H. Morgan, Ph.D., James Moody Ph.D., & Tyler Barret
#3 October 2025

"""
Note: RCall a dependency of the plotting functions used here was removed from the Project toml
      for compilation purposes. If you want to run these tests in Julia, use the command Pkg.add("RCall")
      and load the package.
"""

#	To Do:
#	Create Windows Executable
#	Create R Package that Call Julia Executable

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

#   Setting the DISPLAY
    ENV["DISPLAY"] = ":100"

################
#   PACKAGES   #
################
using CSV
using DataFrames
using EzXML
using ProgressMeter
using Random
using RCall
using StatsBase
using Statistics
using diffusion_sim

######################
#   TEST FUNCTIONS   #
######################

#   Function for Identifying the Mode of a Vector
    function arithmetic_mode(x; na_rm=false)
        """
        arithmetic_mode(x; na_rm=false)
        
        Calculate the mode (most frequent value) of an array.
        
        # Arguments
        - `x`: Input array
        - `na_rm`: If true, remove missing values before calculation (default: false)
        
        # Returns
        - The most frequent value in the array
        """
        
        #   Handle missing values if requested
            if na_rm
                x = filter(!ismissing, x)
            end
        
        #   Get unique values
            vals = unique(x)
        
        #   Count occurrences of each unique value (convert to strings to handle missing)
            x_str = string.(x)
            vals_str = string.(vals)
            counts = [count(==(val), x_str) for val in vals_str]
        
        #   Return the value with maximum count
            return vals[argmax(counts)]
    end

#	Helper Function for plot_parameter_space: assign colors based on infection probability
	function color_assignment(prob_infection::Vector{Float64})
		"""
		Args:
			prob_infection::Vector{Float64}: probability values (0 to 1)
		Returns:
			Vector{String}: color names corresponding to each probability
		Notes:
			Maps probabilities to color gradient from blue to red.
		"""
		
		#	Define color palette
			colorpal = ["slateblue", "royalblue1", "lightblue2", "#FFFFFF",
			           "yellow1", "navajowhite", "lightsalmon", "orange",
			           "darkorange1", "coral", "firebrick1", "brown"]
		
		#	Setting Parameters
			n = length(prob_infection)
			color_vector = Vector{String}(undef, n)
		
		#	Assign colors based on which bin each value falls into
			for i in 1:n
				#	Isolating the Iteration
					val = prob_infection[i]
				
				#	Determine color index based on value ranges
					if val < 0.05
						idx = 1
					elseif val < 0.10
						idx = 2
					elseif val < 0.15
						idx = 3
					elseif val < 0.20
						idx = 4
					elseif val < 0.25
						idx = 5
					elseif val < 0.30
						idx = 6
					elseif val < 0.35
						idx = 7
					elseif val < 0.40
						idx = 8
					elseif val < 0.45
						idx = 9
					elseif val < 0.50
						idx = 10
					elseif val < 0.55
						idx = 11
					else
						idx = 12
					end
				
				#	Populating the Color Vector
					color_vector[i] = colorpal[idx]
			end
		
		#	Return the Color Vector
			return color_vector
	end

#	Plot Parameter Space Using R Graphics
	function plot_parameter_space(results::DataFrame;
	                              edgemax::Int=3,
	                              display::String=":100",
	                              save_path::Union{String,Nothing}=nothing,
	                              show_plot::Bool=true)
		"""
		Args:
			results::DataFrame: output from sim_parm_space function
			edgemax::Int: maximum edge weight (must match sim_parm_space)
			display::String: X11 display target passed to R (e.g., ":100")
			save_path::Union{String,Nothing}: optional file path to save plot
			show_plot::Bool: whether to display plot (default true)
		Returns:
			Nothing
		Notes:
			Creates multi-panel visualization of interaction probabilities.
			Uses R's base graphics via RCall for exact reproduction.
		"""
		
		#	Extract parameters from first row (all rows have same params)
			param_row = results[1, :]
			b_int = param_row.b_int
			b_cxn_symp = param_row.b_cxn_symp
			b_cxn_total = param_row.b_cxn_total
			b_cxn_peers = param_row.b_cxn_peers
			b_close = results[results.edgwgt .== 2, :b_int][1]  # extract b_close from calculation
		
		#	Generate colors for all data
			r_colors = color_assignment(results.prop_inf)
			results.r_colors = r_colors
		
		#	Transfer data to R
			@rput results edgemax b_int b_cxn_symp b_cxn_total b_cxn_peers b_close display
		
		#	Open graphics device if saving
			if save_path !== nothing
				R"""
				pdf($save_path, width=10, height=12)
				"""
			elseif show_plot
				R"""
				#	Set DISPLAY
					Sys.setenv(DISPLAY=display)
					Sys.getenv("DISPLAY")      
					
				#	Creating X11 Window
					X11(type = "cairo", width=10, height=12)
				"""
			end
		
		#	Execute R plotting code
			R"""
			# 	Helper function for plot_parameter_space: Dataplot Graph Ticks
				dataplot_tick_function <- function(major_tick_length=0.035, minor_tick_ratio=0.25){
					# 	Check if Hmisc is Present. If not, install it.
						packages <- c('Hmisc')
						install.packages(setdiff(packages, rownames(installed.packages()))) 
							
					# 	Add Dataplot Style- Through Tick Marks to Plot
						Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = minor_tick_ratio)  
						Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = -minor_tick_ratio)  
						axis(2, tck=1, tck=-major_tick_length, labels = FALSE)
						axis(1, tck=1, tck=-major_tick_length, labels = FALSE)
				}

			# 	Creating Layout
				viz_matrix <- matrix(c(7,7,7,7,7,7,7,7,7,7,
									1,1,1,2,2,2,3,3,3,8,
									1,1,1,2,2,2,3,3,3,8,
									1,1,1,2,2,2,3,3,3,8,
									4,4,4,5,5,5,6,6,6,8,
									4,4,4,5,5,5,6,6,6,8,
									4,4,4,5,5,5,6,6,6,8,
									9,9,9,9,9,9,9,9,9,9), 
									ncol = 10, byrow = TRUE)
				layout(viz_matrix)
			
			# 	Looping Over Visualization Matrix
				for (issymptomatic in seq(0, 1, 1)) {
					# 	Filter for current symptomatic status
						r_panel_data <- results[results$issymptomatic == issymptomatic, ]
					
						for (edgwgt in seq(1, edgemax, 1)) {
							# 	Isolating Cell Data
								r_cell_data <- r_panel_data[r_panel_data$edgwgt == edgwgt, ]
							
							# 	Creating Base Plot
								par(mar=c(3, 4, 2, 2))
								plot(NA, type="n", xlab = "", ylab = "Probability of Interaction", 
									family = 'serif', xlim = c(0, 5), 
									ylim = c(0, 1), las = 1, tck=0.015, xaxt='n', bty='L')

							#	Adding Tick Marks
								axis(1, padj=0.75, tck=0.015) 
              					dataplot_tick_function(0.015, 0.40)

							#	Adding Title
								title(paste('edgwgt =', edgwgt), family='serif', cex.main=1.5, line=-1.15)
								
							# 	Adding Reference Lines
								ref_lines <- seq(0, 1, 0.2)
								for (ref_line in ref_lines) {
									abline(h = ref_line, col = "gray60", lty = 3)
								}
							
							# 	Plotting X-Axis Label
								mtext(side = 1, text = 'Number of Ill Peers', 
									  col = "black", line = 2.55, cex = 0.75, family = 'serif')
							
							# 	Plotting Lines for each ninf level
								ninf_levels <- unique(r_cell_data$ninf)
								for (ninf in ninf_levels) {
									r_line_data <- r_cell_data[r_cell_data$ninf == ninf, ]
									r_line_data <- r_line_data[order(r_line_data$peersinf), ]
									lines(r_line_data$peersinf, r_line_data$prob_act, 
										lty = 1, col = r_line_data$r_colors[1], lwd = 2)
								}
							
							# 	Adding Symptomatic label on rightmost panel
								if (edgwgt == edgemax) {
									mtext(side = 4, text = paste("issymptomatic =", issymptomatic), 
										col = "black", line = 1, cex = 0.75, 
										family = 'serif', font = 2)
								}
						}
				}
			
			#	Adding Plot Title
				par(mar = c(0, 0, 0, 0), bty = 'n')
				plot(0, type = 'n', xlab = ' ', ylab = ' ', 
					cex.axis = 1.3, xaxt = 'n', yaxt = 'n', 
					family = 'serif', las = 1, main = ' ')
				text(x = 1, y = 0.5, 
					labels = "SIR Social Feedback Function Plot", 
					cex = 2.5, family = 'serif', font = 2)
			
			# 	Adding Color Gradient Legend
				unique_props <- sort(unique(results$prop_inf))
				unique_colors <- character(length(unique_props))
				for (i in 1:length(unique_props)) {
					idx <- which(results$prop_inf == unique_props[i])[1]
					unique_colors[i] <- results$r_colors[idx]
				}
			
				legend_image <- as.raster(matrix(rev(unique_colors), ncol = 1))
			
				par(mar = c(0, 0, 0, 0), bty = 'n')
				plot(0, type = 'n', xlab = ' ', ylab = ' ', 
					xaxt = 'n', yaxt = 'n', cex.axis = 1.3, 
					family = 'serif', las = 1, main = ' ')
			
				rasterImage(legend_image, 0.90, -0.75, 1.01, 0.75)
				text(x = 1.15, y = seq(-0.7, 0.7, length.out = 12), 
					labels = sprintf(\"%.2f\", seq(0, 0.55, 0.05)), 
					cex = 1, family = 'serif')
				text(x = 1, y = 0.8, labels = "Prop. Inf.", 
					cex = 1.35, family = 'serif')
			
			# 	Adding Parameter Annotation
				par(mar = c(0, 0, 0, 0), bty = 'n')
				plot(0, type = 'n', xlab = ' ', ylab = ' ', 
					cex.axis = 1.3, xaxt = 'n', yaxt = 'n', 
					family = 'serif', las = 1, main = ' ')
			
				param_text <- paste(
					paste("b_int =", b_int),
					paste("b_cxn_symp =", b_cxn_symp),
					paste("b_cxn_global =", b_cxn_total),
					paste("b_cxn_peers =", b_cxn_peers),
					paste("b_close =", b_close),
					sep = "  "
				)
			
				text(x = 1, y = 0.5, labels = param_text, 
					cex = 1.3, family = 'serif', font = 1)
			"""
		
		#	Close graphics device if saving
			if save_path !== nothing
				R"""
				dev.off()
				"""
				println("Plot saved to: ", save_path)
			end
		
		#	Return nothing
			return nothing
	end
	@doc raw"""
	**Description**
	Creates a multi-panel visualization of SIR model parameter space exploration results.
	Uses R's base graphics via RCall to exactly reproduce the original R plotting code.

	**Usage**
	`plot_parameter_space(results; edgemax=3, display=":100", save_path=nothing, show_plot=true)`

	**Arguments**
	- `results::DataFrame`: Output from sim_parm_space function
	- `edgemax::Int`: Maximum edge weight (must match sim_parm_space call)
	- `display::String`: X11 display target passed to R (e.g., `":100"`)
	- `save_path::Union{String,Nothing}`: Optional file path to save plot as PDF
	- `show_plot::Bool`: Whether to display plot interactively (default true)

	**Details**
	Creates a grid of panels showing:
	- Rows: Symptomatic status (0 = asymptomatic, 1 = symptomatic)
	- Columns: Edge weight values (1 to edgemax)
	- Lines: Different global infection levels (ninf values)
	- X-axis: Number of infected peers
	- Y-axis: Probability of interaction
	- Colors: Proportion of population infected

	The visualization helps assess how social feedback parameters affect contact patterns
	under different epidemic and network conditions.

	**Value**
	Returns nothing. Creates plot as side effect.

	**Examples**
	```julia
	# Generate parameter space
	results = sim_parm_space(-0.1, -1.5, -3.5, -0.5, 1.0, -0.1, 3, 5, 500)
	
	# Display plot interactively
	plot_parameter_space(results, edgemax=3, display=":100")
	
	# Save to file
	plot_parameter_space(results, edgemax=3, save_path="sir_parameters.pdf", show_plot=false)
	```

	**See Also**
	`sim_parm_space` for generating the data to visualize
	""" plot_parameter_space


#	Parameter Space Exploration for SIR Diffusion Simulation
	function sim_parm_space(b_int::Float64, b_cxn_symp::Float64, b_cxn_global::Float64, 
	                        b_cxn_peers::Float64, b_close::Float64, b_cls_x_smp::Float64,
	                        edgemax::Int, maxdeg::Int, n::Int)
		"""
		Args:
			b_int::Float64: baseline activity/interaction parameter
			b_cxn_symp::Float64: activity reduction when symptomatic
			b_cxn_global::Float64: activity reduction from global infection level
			b_cxn_peers::Float64: activity reduction from infected peers
			b_close::Float64: activity boost from edge weight (connection strength)
			b_cls_x_smp::Float64: interaction term between edge weight and symptomatic status
			edgemax::Int: maximum edge weight in network
			maxdeg::Int: maximum degree (number of connections) a node can have
			n::Int: population size for scaling global infection impact
		Returns:
			DataFrame: containing prob_act, parameter values, and conditions for each combination
		Notes:
			Explores parameter space for social feedback in epidemic spread.
			Uses logistic transformation to calculate interaction probabilities.
			Fixed population of 500 for prop_inf calculation.
		"""
		
		#	Calculate total number of parameter combinations
			n_symptomatic = 2  # 0 and 1
			n_inf_levels = length(1:50:300)  # 6 levels
			n_rows = n_symptomatic * n_inf_levels * edgemax * (maxdeg + 1)
		
		#	Pre-allocate arrays for DataFrame
			prob_act = Vector{Float64}(undef, n_rows)
			issymptomatic = Vector{Int}(undef, n_rows)
			ninf = Vector{Int}(undef, n_rows)
			edgwgt = Vector{Int}(undef, n_rows)
			peersinf = Vector{Int}(undef, n_rows)
			b_int_vec = Vector{Float64}(undef, n_rows)
			b_cxn_symp_vec = Vector{Float64}(undef, n_rows)
			b_cxn_total_vec = Vector{Float64}(undef, n_rows)
			b_cxn_peers_vec = Vector{Float64}(undef, n_rows)
		
		#	Parameter aliasing for clarity
			b_cxn_total = b_cxn_global  # consistent with SAS naming
			
		#	Iterate over parameter combinations with index tracking
			idx = 1
			for issymp in 0:1
				for n_infected in 1:50:300  # 1, 51, 101, 151, 201, 251
					for edge_weight in 1:edgemax
						for peers_infected in 0:maxdeg
							
							#	Calculate log-odds of activity
								lwact = b_int + 
								       (issymp * b_cxn_symp) + 
								       (b_close * edge_weight) + 
								       (b_cxn_peers * peers_infected) + 
								       b_cxn_total * (n_infected / (n/3)) + 
								       b_cls_x_smp * edge_weight * issymp
							
							#	Transform to probability and populate arrays
								prob_act[idx] = exp(lwact) / (1 + exp(lwact))
								issymptomatic[idx] = issymp
								ninf[idx] = n_infected
								edgwgt[idx] = edge_weight
								peersinf[idx] = peers_infected
								b_int_vec[idx] = b_int
								b_cxn_symp_vec[idx] = b_cxn_symp
								b_cxn_total_vec[idx] = b_cxn_total
								b_cxn_peers_vec[idx] = b_cxn_peers
								
							#	Increment index
								idx += 1
						end
					end
				end
			end
		
		#	Construct DataFrame from pre-allocated arrays
			results = DataFrame(
				prob_act = prob_act,
				issymptomatic = issymptomatic,
				ninf = ninf,
				edgwgt = edgwgt,
				peersinf = peersinf,
				b_int = b_int_vec,
				b_cxn_symp = b_cxn_symp_vec,
				b_cxn_total = b_cxn_total_vec,
				b_cxn_peers = b_cxn_peers_vec
			)
		
		#	Add proportion infected column
			results.prop_inf = results.ninf ./ 500
			
		#	Return completed dataframe
			return results
	end
	@doc raw"""
	**Description**
	Explores the parameter space for an SIR epidemic model with social feedback mechanisms. 
	Calculates interaction probabilities based on individual symptomatic status, local peer 
	infection, global infection levels, and network connection strengths.

	**Usage**
	`sim_parm_space(b_int, b_cxn_symp, b_cxn_global, b_cxn_peers, b_close, b_cls_x_smp, edgemax, maxdeg, n)`

	**Arguments**
	- `b_int::Float64`: Baseline interaction/activity parameter (intercept in logistic model)
	- `b_cxn_symp::Float64`: Reduction in activity when individual is symptomatic
	- `b_cxn_global::Float64`: Activity reduction based on global infection prevalence  
	- `b_cxn_peers::Float64`: Activity reduction based on number of infected peers
	- `b_close::Float64`: Activity boost from edge weight (strength of connection)
	- `b_cls_x_smp::Float64`: Interaction effect between edge weight and symptomatic status
	- `edgemax::Int`: Maximum edge weight value to explore
	- `maxdeg::Int`: Maximum node degree (number of peers) to explore
	- `n::Int`: Population size for scaling global infection impact (typically 500)

	**Details**
	The function systematically explores combinations of:
	- Symptomatic status (0 or 1)
	- Number infected in population (1 to 300 by steps of 50)
	- Edge weights (1 to edgemax)
	- Number of infected peers (0 to maxdeg)

	For each combination, calculates the log-odds of interaction:
	```
	lwact = b_int + (issymptomatic * b_cxn_symp) + (b_close * edgwgt) + 
	        (b_cxn_peers * peersinf) + b_cxn_total * (ninf/(n/3)) + 
	        b_cls_x_smp * edgwgt * issymptomatic
	```
	
	Then transforms to probability using the logistic function: `p = exp(lwact)/(1 + exp(lwact))`

	The term `ninf/(n/3)` scales the global infection effect, assuming roughly 1/3 of 
	population could be infected at peak.

	**Value**
	Returns a `DataFrame` with columns:
	- `prob_act`: Probability of interaction/activity
	- `issymptomatic`: Whether individual is symptomatic (0/1)
	- `ninf`: Number infected in population
	- `edgwgt`: Edge weight value
	- `peersinf`: Number of infected peers
	- `b_int`, `b_cxn_symp`, `b_cxn_total`, `b_cxn_peers`: Parameter values used
	- `prop_inf`: Proportion infected (ninf/500)

	**Examples**
	```julia
	# Explore parameter space with moderate social feedback
	results = sim_parm_space(-0.1, -1.5, -3.5, -0.5, 1.0, -0.1, 3, 5, 500)
	
	# Filter for specific conditions
	using DataFramesMeta
	symptomatic_results = @subset(results, :issymptomatic .== 1)
	```

	**See Also**
	Related epidemic simulation functions that use these parameters for agent-based modeling.
	""" sim_parm_space

#	Test Suite for sim_parm_space Function
	function test_sim_parm_space(; verbose::Bool=true, export_csv::Bool=false)
		"""
		Args:
			verbose::Bool: whether to print detailed output (default true)
			export_csv::Bool: whether to export CSV files (default false)
		Returns:
			NamedTuple: test results and generated data
		Notes:
			Comprehensive test suite for parameter space exploration.
			Includes validation, sensitivity analysis, and export options.
		"""
		
		#	Initialize test results storage
			test_results = Dict{String, Any}()
			
		#	Run basic parameter exploration
			if verbose
				println("="^60)
				println("Running sim_parm_space tests")
				println("="^60)
			end
			
			results = sim_parm_space(
				-0.1,   # b_int: baseline activity
				-1.5,   # b_cxn_symp: activity reduction when symptomatic
				-3.5,   # b_cxn_global: global infection impact
				-0.5,   # b_cxn_peers: peer infection impact
				1.0,    # b_close: edge weight boost
				-0.1,   # b_cls_x_smp: interaction term
				3,      # edgemax: max edge weight
				5,      # maxdeg: max node degree
				500     # n: population size
			)
			
			test_results["main_results"] = results
			test_results["n_rows"] = nrow(results)
			
		#	Test 1: Verify dimensions
			expected_rows = 2 * 6 * 3 * 6  # symptomatic * ninf_levels * edgemax * (maxdeg+1)
			dimension_test = nrow(results) == expected_rows
			test_results["dimension_test"] = dimension_test
			
			if verbose
				println("\nTest 1: Dimension Check")
				println("  Expected rows: ", expected_rows)
				println("  Actual rows: ", nrow(results))
				println("  ✓ Pass: ", dimension_test)
			end
			
		#	Test 2: Probability bounds
			prob_bounds_test = all(0 .≤ results.prob_act .≤ 1)
			test_results["prob_bounds_test"] = prob_bounds_test
			
			if verbose
				println("\nTest 2: Probability Bounds")
				println("  Min probability: ", round(minimum(results.prob_act), digits=6))
				println("  Max probability: ", round(maximum(results.prob_act), digits=6))
				println("  ✓ Pass: ", prob_bounds_test)
			end
			
		#	Test 3: Manual calculation verification
			test_row = results[100, :]
			lwact = test_row.b_int + 
			       (test_row.issymptomatic * test_row.b_cxn_symp) + 
			       (1.0 * test_row.edgwgt) + 
			       (test_row.b_cxn_peers * test_row.peersinf) + 
			       test_row.b_cxn_total * (test_row.ninf / (500/3)) + 
			       (-0.1) * test_row.edgwgt * test_row.issymptomatic
			
			manual_prob = exp(lwact) / (1 + exp(lwact))
			calculation_test = isapprox(test_row.prob_act, manual_prob, rtol=1e-10)
			test_results["calculation_test"] = calculation_test
			
			if verbose
				println("\nTest 3: Calculation Verification")
				println("  Row 100 stored value: ", round(test_row.prob_act, digits=6))
				println("  Manual calculation:   ", round(manual_prob, digits=6))
				println("  ✓ Pass: ", calculation_test)
			end
			
		#	Test 4: Monotonicity checks
			#	Check: increasing peersinf should decrease prob_act (given negative b_cxn_peers)
			mono_subset = filter(row -> 
				row.issymptomatic == 0 && 
				row.ninf == 51 && 
				row.edgwgt == 2, 
				results)
			
			mono_test = issorted(mono_subset.prob_act, rev=true)  # Should be decreasing
			test_results["monotonicity_test"] = mono_test
			
			if verbose
				println("\nTest 4: Monotonicity Check")
				println("  Testing: prob_act decreases as peersinf increases")
				println("  (fixed: issymptomatic=0, ninf=51, edgwgt=2)")
				println("  ✓ Pass: ", mono_test)
			end
			
		#	Test 5: Summary statistics
			grouped = groupby(results, :issymptomatic)
			summary_stats = combine(grouped,
				:prob_act => mean => :mean_prob,
				:prob_act => std => :std_prob,
				:prob_act => minimum => :min_prob,
				:prob_act => maximum => :max_prob
			)
			test_results["summary_stats"] = summary_stats
			
			if verbose
				println("\nTest 5: Summary Statistics by Symptomatic Status")
				for row in eachrow(summary_stats)
					println("  Symptomatic = ", row.issymptomatic, ":")
					println("    Mean: ", round(row.mean_prob, digits=4), 
						   " (SD: ", round(row.std_prob, digits=4), ")")
					println("    Range: [", round(row.min_prob, digits=4), 
						   ", ", round(row.max_prob, digits=4), "]")
				end
			end
			
		#	Test 6: Extreme scenarios
			baseline = filter(row -> 
				row.issymptomatic == 0 && 
				row.ninf == 1 && 
				row.edgwgt == 1 && 
				row.peersinf == 0, 
				results)[1, :]
			
			worst_case = filter(row -> 
				row.issymptomatic == 1 && 
				row.ninf == 251 && 
				row.edgwgt == 3 && 
				row.peersinf == 5, 
				results)[1, :]
			
			reduction_factor = worst_case.prob_act / baseline.prob_act
			test_results["baseline_prob"] = baseline.prob_act
			test_results["worst_case_prob"] = worst_case.prob_act
			test_results["reduction_factor"] = reduction_factor
			
			if verbose
				println("\nTest 6: Extreme Scenarios")
				println("  Baseline (minimal risk):  ", round(baseline.prob_act, digits=4))
				println("  Worst case (maximum risk): ", round(worst_case.prob_act, digits=4))
				println("  Reduction factor: ", round(reduction_factor, digits=4))
			end
			
		#	Export if requested
			if export_csv
				CSV.write("sir_parameter_space.csv", results)
				test_results["export_full"] = true
				
				#	Create visualization subset
				viz_subset = filter(row -> 
					row.ninf ∈ [1, 101, 201] && 
					row.edgwgt ∈ [1, 3], 
					results)
				CSV.write("sir_parameter_space_viz.csv", viz_subset)
				test_results["export_viz"] = true
				
				if verbose
					println("\nExport Complete:")
					println("  Full dataset: sir_parameter_space.csv (", nrow(results), " rows)")
					println("  Viz subset: sir_parameter_space_viz.csv (", nrow(viz_subset), " rows)")
				end
			end
			
		#	Summary
			all_tests_passed = dimension_test && prob_bounds_test && 
			                  calculation_test && mono_test
			test_results["all_passed"] = all_tests_passed
			
			if verbose
				println("\n" * "="^60)
				println("Test Summary: ", all_tests_passed ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗")
				println("="^60)
			end
			
		#	Return test results as NamedTuple
			return (
				results = results,
				passed = all_tests_passed,
				test_details = test_results,
				summary_stats = summary_stats
			)
	end
	@doc raw"""
	**Description**
	Comprehensive test suite for the sim_parm_space function. Validates calculations,
	checks probability bounds, verifies monotonicity, and optionally exports results.

	**Usage**
	`test_sim_parm_space(; verbose=true, export_csv=false)`

	**Arguments**
	- `verbose::Bool`: Print detailed test output (default true)
	- `export_csv::Bool`: Export results to CSV files (default false)

	**Details**
	Runs the following tests:
	1. Dimension verification (expected vs actual row count)
	2. Probability bounds (all values in [0,1])
	3. Manual calculation verification
	4. Monotonicity checks (decreasing probability with more infected peers)
	5. Summary statistics by symptomatic status
	6. Extreme scenario comparison

	**Value**
	Returns a NamedTuple containing:
	- `results`: Full DataFrame from sim_parm_space
	- `passed`: Boolean indicating if all tests passed
	- `test_details`: Dictionary with detailed test results
	- `summary_stats`: Summary statistics by symptomatic status

	**Examples**
	```julia
	# Run tests with output
	test_results = test_sim_parm_space(verbose=true, export_csv=false)
	
	# Silent testing for automated checks
	test_results = test_sim_parm_space(verbose=false, export_csv=false)
	@assert test_results.passed "Tests failed!"
	
	# Generate CSV files for R analysis
	test_results = test_sim_parm_space(verbose=true, export_csv=true)
	```
	""" test_sim_parm_space

#	Quick test function for specific parameter combinations
	function quick_param_test(b_int::Float64, b_cxn_symp::Float64)
		"""
		Args:
			b_int::Float64: baseline activity parameter to test
			b_cxn_symp::Float64: symptomatic reduction parameter to test
		Returns:
			DataFrame: subset of results for quick inspection
		Notes:
			Runs sim_parm_space with reduced parameter space for quick testing.
		"""
		
		#	Run with reduced parameter space
			results = sim_parm_space(
				b_int,
				b_cxn_symp,
				-2.0,   # moderate global effect
				-0.3,   # moderate peer effect
				0.5,    # moderate edge boost
				-0.05,  # small interaction term
				2,      # just 2 edge weights
				3,      # just 3 max peers
				500
			)
			
		#	Return summary
			println("Quick test with b_int=", b_int, ", b_cxn_symp=", b_cxn_symp)
			println("Generated ", nrow(results), " rows")
			println("Probability range: [", 
				round(minimum(results.prob_act), digits=4), ", ",
				round(maximum(results.prob_act), digits=4), "]")
			
			return results
	end

#	Test Function for SIR Model
	function test_sir_model(network_data::NamedTuple;
	                        seed::Int=42,
	                        verbose::Bool=true)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
			seed::Int: random seed for reproducibility (used for tie-breakers; seeds fixed below)
			verbose::Bool: print detailed output
		Returns:
			Dict: test results and metrics
		Notes:
			R-aligned test: 5 fixed seeds, inf_r=0.33, rec_t=14, maxtime=200, p_symp=0.5,
			b_int=-0.1, b_close=1.0, b_cxn_peers=-0.5, b_cxn_total=-3.5, b_cxn_symp=-1.5, b_cls_x_smp=-0.1.
		"""

		#	Set random seed (for reproducibility of any RNG inside sirdif)
			Random.seed!(seed)

		#	Use matrix format directly from sim_prep (sirdif will convert internally)
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)	# Get n from matrix dimensions

		#	Set simulation parameters (R trial settings)
			params = Dict(
				"inf_r"        => 0.33,		# infection rate
				"rec_t"        => 14,			# recovery time
				"maxtime"      => 200,		# simulation duration
				"p_symp"       => 0.5,		# probability symptomatic
				"b_int"        => -0.1,		# baseline
				"b_close"      => 1.0,		# edge weight effect
				"b_cxn_peers"  => -0.5,		# peer effect
				"b_cxn_total"  => -3.5,		# global effect
				"b_cxn_symp"   => -1.5,		# symptomatic effect
				"b_cls_x_smp"  => -0.1		# interaction term
			)

		#	Set fixed seed nodes (R trial): ensure within bounds 1:n
			raw_seeds = [5, 10, 15, 20, 25]
			seed_nodes = [s for s in raw_seeds if 1 ≤ s ≤ n]
			if isempty(seed_nodes)
				error("Provided seed nodes $(raw_seeds) are out of bounds for n = $n.")
			end

			if verbose
				println("="^60)
				println("SIR Model Test (R Settings)")
				println("="^60)
				println("Network size: $n nodes")
				println("Seed nodes: $(seed_nodes)")
				println("Parameters:")
				for (k, v) in params
					println("  $k: $v")
				end
			end

		#	Run simulation (SIR: permanent immunity)
			results = sirdif(
				alst, vlst, seed_nodes,
				params["inf_r"], params["rec_t"], params["maxtime"], params["p_symp"],
				params["b_int"], params["b_close"], params["b_cxn_peers"],
				params["b_cxn_total"], params["b_cxn_symp"], params["b_cls_x_smp"];
				transmission_method = :weighted,
				immunity_duration   = nothing
			)

		#	Extract results
			infection_log = results["infection_log"]
			runtime = results["total_time"]

		#	Calculate key metrics
			final_row = infection_log[end, :]
			final_size = final_row[3]							# prop_ever_infected
			peak_infected = maximum(infection_log[:, 2])		# max n_infected
			peak_time = findfirst(x -> x == peak_infected, infection_log[:, 2])

		#	Validation checks
			checks = Dict{String, Bool}()

			#	Check 1: Monotonic increase in ever infected
				ever_infected = infection_log[:, 3]
				checks["monotonic_ever_infected"] = all(diff(ever_infected) .>= 0)

			#	Check 2: Conservation (infected + recovered ≤ n)
				total_affected = infection_log[:, 2] .+ infection_log[:, 5]
				checks["conservation"] = all(total_affected .<= n)

			#	Check 3: No resurrection (recovered stay recovered in SIR)
				recovered = infection_log[:, 5]
				checks["no_resurrection"] = all(diff(recovered) .>= 0)

			#	Check 4: Epidemic ends (no infected at end or before maxtime)
				checks["epidemic_ends"] = infection_log[end, 2] == 0

		#	Summary statistics
			if verbose
				println("\n" * "="^40)
				println("Results:")
				println("  Runtime: $(round(runtime, digits=3)) seconds")
				println("  Final size: $(round(final_size * 100, digits=2))%  (≈ $(round(final_size * n; digits=0)) of $n)")
				println("  Peak infected: $peak_infected at time $peak_time")
				println("  Epidemic duration: $(size(infection_log, 1) - 1) days")

				println("\nValidation Checks:")
				for (check, passed) in checks
					status = passed ? "✓ PASS" : "✗ FAIL"
					println("  $check: $status")
				end

				all_passed = all(values(checks))
				println("\nOverall: " * (all_passed ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗"))
			end

		#	Return comprehensive results
			return Dict(
				"results" => results,
				"metrics" => Dict(
					"final_size" => final_size,
					"peak_infected" => peak_infected,
					"peak_time" => peak_time,
					"duration" => size(infection_log, 1) - 1
				),
				"checks" => checks,
				"params" => params,
				"seeds"  => seed_nodes
			)
	end

#	Test Function for SIRS Model
	function test_sirs_model(network_data::NamedTuple;
	                         immunity_days::Int=30,
	                         seed::Int=42,
	                         verbose::Bool=true)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep (matrix form: .alst, .vlst)
			immunity_days::Int: duration of immunity after recovery
			seed::Int: random seed for reproducibility
			verbose::Bool: print detailed output
		Returns:
			Dict: test results and metrics
		Notes:
			R-aligned SIRS test:
			- Seeds = [5, 10, 15, 20, 25]
			- inf_r=0.33, rec_t=14, maxtime=400, p_symp=0.5
			- b_int=-0.1, b_close=1.0, b_cxn_peers=-0.25, b_cxn_total=-2.0,
			  b_cxn_symp=-1.5, b_cls_x_smp=-0.1
			- SIRS enabled via immunity_duration=immunity_days

			Infection log columns expected from `sirdif`:
			[ time, NI, prop_ever_inf, prop_cur_inf, NR, prop_recovered, NImm ]
		"""

		#	Set random seed
			Random.seed!(seed)

		#	Use matrix format directly from sim_prep (sirdif will convert internally)
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)

		#	Set simulation parameters (R trial settings)
			params = Dict(
				"inf_r"        => 0.33,
				"rec_t"        => 14,
				"maxtime"      => 400,
				"p_symp"       => 0.5,
				"b_int"        => -0.1,
				"b_close"      => 1.0,
				"b_cxn_peers"  => -0.25,
				"b_cxn_total"  => -2.0,
				"b_cxn_symp"   => -1.5,
				"b_cls_x_smp"  => -0.1
			)

		#	Set fixed seed nodes (R trial): ensure within bounds 1:n
			raw_seeds = [5, 10, 15, 20, 25]
			seed_nodes = [s for s in raw_seeds if 1 ≤ s ≤ n]
			if isempty(seed_nodes)
				error("Provided seed nodes $(raw_seeds) are out of bounds for n = $n.")
			end

			if verbose
				println("="^60)
				println("SIRS Model Test (R Settings)")
				println("="^60)
				println("Network size: $n nodes")
				println("Seed nodes: $(seed_nodes)")
				println("Immunity duration: $immunity_days days")
				println("Parameters:")
				for (k, v) in params
					println("  $k: $v")
				end
			end

		#	Run simulation (SIRS via immunity_duration)
			results = sirdif(
				alst, vlst, seed_nodes,
				params["inf_r"], params["rec_t"], params["maxtime"], params["p_symp"],
				params["b_int"], params["b_close"], params["b_cxn_peers"],
				params["b_cxn_total"], params["b_cxn_symp"], params["b_cls_x_smp"];
				transmission_method = :weighted,
				immunity_duration   = immunity_days
			)

		#	Extract results
			infection_log = results["infection_log"]
			runtime = results["total_time"]

		#	Calculate key metrics
			final_row = infection_log[end, :]
			final_size = final_row[3]						# cumulative proportion ever infected
			peak_infected = maximum(infection_log[:, 2])
			peak_time = findfirst(x -> x == peak_infected, infection_log[:, 2])

		#	Current immune count series (from sirdif column 7)
			NImm_series = infection_log[:, 7]
			peak_immune = maximum(NImm_series)
			peak_immune_time = findfirst(x -> x == peak_immune, NImm_series)

		#	Look for reinfection waves (robust peak detection)
			infected_series = infection_log[:, 2]

			#	Smooth lightly (5-day moving average) to suppress single-day noise
				smoothed = similar(infected_series)
				window = 5
				halfwin = fld(window, 2)
				for i in eachindex(infected_series)
					l = max(1, i - halfwin)
					r = min(length(infected_series), i + halfwin)
					smoothed[i] = sum(infected_series[l:r]) / (r - l + 1)
				end

			#	Peak selection rules
				min_gap = immunity_days + params["rec_t"]				# required separation between waves
				peak_height_thresh = max(5.0, 0.05 * maximum(smoothed))	# ≥5 cases or ≥5% of peak

			#	Find candidate local maxima (on smoothed series)
				candidates = Int[]
				for i in 2:(length(smoothed)-1)
					if (smoothed[i] > smoothed[i-1]) && (smoothed[i] > smoothed[i+1]) && (smoothed[i] ≥ peak_height_thresh)
						push!(candidates, i)
					end
				end

			#	Enforce minimum gap: keep taller peak if two are too close
				peaks = Int[]
				for idx in candidates
					if isempty(peaks) || (idx - peaks[end] ≥ min_gap)
						push!(peaks, idx)
					else
						if smoothed[idx] > smoothed[peaks[end]]
							peaks[end] = idx
						end
					end
				end
				n_waves = length(peaks)

		#	Validation checks
			checks = Dict{String, Bool}()

			#	Check 1: State bounds (SIRS-aware)
				checks["infected_within_bounds"] = all((infection_log[:, 2] .>= 0) .& (infection_log[:, 2] .<= n))
				checks["immune_within_bounds"]   = all((NImm_series .>= 0) .& (NImm_series .<= n))
				checks["sirS_conservation"]      = all((infection_log[:, 2] .+ NImm_series) .<= n)

			#	Check 2: Reinfection timing (using filtered peaks)
				if n_waves > 1
					wave_gaps = diff(peaks)
					checks["reinfection_timing"] = all(wave_gaps .>= min_gap)
				else
					checks["reinfection_timing"] = true
				end

			#	Check 3: Reasonable total (cumulative measure under SIRS)
				max_reinfections = params["maxtime"] ÷ (params["rec_t"] + immunity_days)
				checks["reasonable_total"] = final_size <= max_reinfections

		#	Summary statistics
			if verbose
				println("\n" * "="^40)
				println("Results:")
				println("  Runtime: $(round(runtime, digits=3)) seconds")
				println("  Final size: $(round(final_size * 100, digits=2))%  (cumulative)")
				println("  Peak infected: $peak_infected at time $peak_time")
				println("  Peak immune: $peak_immune at time $peak_immune_time")
				println("  Number of waves: $n_waves")
				println("  Epidemic duration: $(size(infection_log, 1) - 1) days")

				println("\nValidation Checks:")
				for (check, passed) in checks
					status = passed ? "✓ PASS" : "✗ FAIL"
					println("  $check: $status")
				end

				all_passed = all(values(checks))
				println("\nOverall: " * (all_passed ? "ALL TESTS PASSED ✓" : "SOME TESTS FAILED ✗"))
			end

		#	Return comprehensive results
			return Dict(
				"results" => results,
				"metrics" => Dict(
					"final_size" => final_size,
					"peak_infected" => peak_infected,
					"peak_time" => peak_time,
					"peak_immune" => peak_immune,
					"peak_immune_time" => peak_immune_time,
					"duration" => size(infection_log, 1) - 1,
					"n_waves" => n_waves
				),
				"checks" => checks,
				"params" => params,
				"seeds"  => seed_nodes
			)
	end

#	Helper Function for replicate_sas_simulation: build design matrix of coefficients
	function replicate_sas_simulation_build_design(bproot::Float64, bgroot::Float64, bsroot::Float64)
		"""
		Args:
			bproot::Float64: root value for b_cxn_peer  (SAS bproot)
			bgroot::Float64: root value for b_cxn_global (SAS bgroot)
			bsroot::Float64: root value for b_cxn_symp  (SAS bsroot)
		Returns:
			Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}:
				(b_peer_v, b_glob_v, b_self_v) each length 96 (4×4×6)
		Notes:
			Mirrors the SAS nested loops:
			  b_peer ∈ {0, bproot, bproot-0.25, bproot-0.50}
			  b_glob ∈ {0, bgroot, bgroot-1.5, bgroot-3.0}
			  b_self ∈ {0, bsroot, bsroot-0.5, bsroot-1.0, bsroot-1.5, bsroot-2.0}
			Order is deterministic: peer major, then global, then self.
		"""

		#	Cardinalities and allocation
			n_peer, n_glob, n_self = 4, 4, 6
			total_coef = n_peer * n_glob * n_self
			b_peer_v  = Vector{Float64}(undef, total_coef)
			b_glob_v  = Vector{Float64}(undef, total_coef)
			b_self_v  = Vector{Float64}(undef, total_coef)

		#	Fill sequences
			r = 0
			for b_peer in 1:n_peer
				b_cxn_peer = (b_peer == 1) ? 0.0 :
				              (b_peer == 2) ? bproot :
				              bproot - 0.25 * (b_peer - 2)

				for b_glob in 1:n_glob
					b_cxn_global = (b_glob == 1) ? 0.0 :
					               (b_glob == 2) ? bgroot :
					               bgroot - 1.5 * (b_glob - 2)

					for b_s in 1:n_self
						b_self = (b_s == 1) ? 0.0 :
						         (b_s == 2) ? bsroot :
						         bsroot - 0.5 * (b_s - 2)

						r += 1
						b_peer_v[r] = b_cxn_peer
						b_glob_v[r] = b_cxn_global
						b_self_v[r] = b_self
					end
				end
			end

		#	Return design vectors
			return (b_peer_v, b_glob_v, b_self_v)
	end
		
#	Helper Function for compare_parameter_sweeps: replicate_sas_parameter_sweep
	function replicate_sas_simulation(
		alst, vlst,           # fixed network (matrix or vector-of-vectors; must match)
		inf_r::Float64,       # base transmission prob (SAS bt)
		rec_t::Int,           # recovery time (days)
		maxtime::Int,         # horizon
		num_iter::Int         # iterations per design row
	)
		"""
		Args:
			alst: adjacency (matrix or vector-of-vectors)
			vlst: edge weights (matching `alst`)
			inf_r::Float64: base transmission probability (SAS `bt`)
			rec_t::Int: recovery time in days
			maxtime::Int: maximum simulation days (≥ 1)
			num_iter::Int: iterations per coefficient setting (≥ 1)
		Returns:
			DataFrame: long-form results with **all varied parameters recorded per row**
		Notes:
			- Builds a coefficient **design matrix**, then loops rows → iterations.
			- Records: `b_cxn_peer`, `b_cxn_global`, `b_self` (varied), plus constants
			(`b_int`, `b_close`, `b_cls_x_smp`, `p_symp`) and run ids (`coef_row`, `itter`, `seednode`).
			- Preallocates to `96 * num_iter * (maxtime+1)` rows (upper bound), then trims.
			- Progress bar ticks once per (design row × iteration).
		"""

		#	Validation
			if maxtime < 1
				throw(ArgumentError("maxtime must be ≥ 1"))
			end
			if num_iter < 1
				throw(ArgumentError("num_iter must be ≥ 1"))
			end

		#	Held-constant parameters (SAS defaults)
			params = Dict{Symbol,Any}(
				:b_int             => -0.1,
				:b_close           =>  1.0,
				:b_cxn_peers_root  => -0.5,
				:b_cxn_global_root => -3.5,
				:b_cxn_symp_root   => -1.5,
				:b_cls_x_smp       => -0.1,
				:p_symp            =>  0.75,
				:progress_desc     => "SAS IML-style simulation (design-driven)",
			)

		#	Build coefficient design (4×4×6 = 96 rows)
			b_peer_v, b_glob_v, b_self_v = replicate_sas_simulation_build_design(
				params[:b_cxn_peers_root], params[:b_cxn_global_root], params[:b_cxn_symp_root],
			)
			total_rows = length(b_peer_v)  # 96

		#	Progress bar
			p = ProgressMeter.Progress(total_rows * num_iter; desc = params[:progress_desc], dt = 0.1)

		#	Upper bound on output rows (include t=0)
			max_out = total_rows * num_iter * (maxtime + 1)

		#	Preallocate outputs (record all varied params per record)
			coef_row_v     = Vector{Int}(undef, max_out)
			itter_v        = Vector{Int}(undef, max_out)
			seednode_v     = Vector{Int}(undef, max_out)

		#	coefficients (varied) + constants (for traceability)
			b_cxn_peer_v   = Vector{Float64}(undef, max_out)
			b_cxn_global_v = Vector{Float64}(undef, max_out)
			b_self_v_out   = Vector{Float64}(undef, max_out)
			b_int_v        = Vector{Float64}(undef, max_out)
			b_close_v      = Vector{Float64}(undef, max_out)
			b_cls_x_smp_v  = Vector{Float64}(undef, max_out)
			p_symp_v       = Vector{Float64}(undef, max_out)
			beta_v         = Vector{Float64}(undef, max_out)

		#	outcomes (per time step)
			time_v         = Vector{Float64}(undef, max_out)
			propeverinf_v  = Vector{Float64}(undef, max_out)
			propcurrinf_v  = Vector{Float64}(undef, max_out)
			proprec_v      = Vector{Float64}(undef, max_out)

		#	Write cursor and node count
			w = 0
			n_nodes = isa(alst, Matrix) ? size(alst, 1) : length(alst)

		#	Main loops: design rows → iterations
			for i in 1:total_rows
				b_peers = b_peer_v[i]
				b_glob  = b_glob_v[i]
				b_self  = b_self_v[i]

				for it in 1:num_iter
					#	Random single seed (match SAS sampling of one node)
						seednode = rand(1:n_nodes)
						infectedp = [seednode]

					#	Run sirdif with this row's coefficients
						result = sirdif(
							alst, vlst, infectedp,
							inf_r, rec_t, maxtime,
							params[:p_symp],
							params[:b_int], params[:b_close],
							b_peers, b_glob, b_self, params[:b_cls_x_smp];
							transmission_method = :weighted,
							immunity_duration   = nothing,
						)

					#	Copy trajectory rows (now a DataFrame)
						logdf = result["infection_log"]::DataFrame
						T     = nrow(logdf)

						@inbounds for k in 1:T
							w += 1
							coef_row_v[w]     = i
							itter_v[w]        = it
							seednode_v[w]     = seednode

							b_cxn_peer_v[w]   = b_peers
							b_cxn_global_v[w] = b_glob
							b_self_v_out[w]   = b_self
							b_int_v[w]        = params[:b_int]
							b_close_v[w]      = params[:b_close]
							b_cls_x_smp_v[w]  = params[:b_cls_x_smp]
							p_symp_v[w]       = params[:p_symp]
							beta_v[w]         = inf_r

							time_v[w]         = logdf.time[k]
							propeverinf_v[w]  = logdf.prop_ever_infected[k]
							propcurrinf_v[w]  = logdf.prop_currently_infected[k]
							proprec_v[w]      = logdf.prop_recovered[k]
						end

					ProgressMeter.next!(p)
				end
			end

		#	Trim and assemble DataFrame
			last = w
			return DataFrame(
				coef_row     = coef_row_v[1:last],
				itter        = itter_v[1:last],
				seednode     = seednode_v[1:last],

				b_cxn_peer   = b_cxn_peer_v[1:last],
				b_cxn_global = b_cxn_global_v[1:last],
				b_self       = b_self_v_out[1:last],
				b_int        = b_int_v[1:last],
				b_close      = b_close_v[1:last],
				b_cls_x_smp  = b_cls_x_smp_v[1:last],
				p_symp       = p_symp_v[1:last],
				beta         = beta_v[1:last],

				time         = time_v[1:last],
				propeverinf  = propeverinf_v[1:last],
				propcurrinf  = propcurrinf_v[1:last],
				proprec      = proprec_v[1:last],
			)
	end
	@doc raw"""
	**Description**
	Replicates the SAS PROC IML simulation on a fixed input network by (1) building a
	coefficient design matrix and (2) iterating over design rows and iterations. For each
	run, one random seed node is chosen and `sirdif` is executed; all **varied parameters**
	are recorded on every output row to enable aggregation (means/medians) and plotting.

	**Usage**
	`replicate_sas_simulation(alst, vlst, inf_r, rec_t, maxtime, num_iter)`

	**Arguments**
	- `alst`, `vlst`: Network adjacency and weights (matrix or vector-of-vectors; same form).
	- `inf_r::Float64`: Base transmission probability per contact.
	- `rec_t::Int`: Recovery time in days.
	- `maxtime::Int`: Maximum simulation days (≥ 1).
	- `num_iter::Int`: Iterations per coefficient setting (≥ 1).

	**Details**
	- Internal parameters (held constant): `b_int=-0.1`, `b_close=1.0`, `b_cls_x_smp=-0.1`, `p_symp=0.75`.
	- Coefficient design matrix (96 rows) mirrors SAS sequences:
	- `b_cxn_peer ∈ {0, bproot, bproot-0.25, bproot-0.50}`, with `bproot=-0.5`
	- `b_cxn_global ∈ {0, bgroot, bgroot-1.5, bgroot-3.0}`, with `bgroot=-3.5`
	- `b_self ∈ {0, bsroot, bsroot-0.5, bsroot-1.0, bsroot-1.5, bsroot-2.0}`, with `bsroot=-1.5`
	- For each `(design row × iteration)` the full `sirdif` trajectory is copied to the
	result. A progress bar advances once per run.

	**Value**
	Returns a long-form `DataFrame` with columns:
	- Run identifiers: `coef_row::Int`, `itter::Int`, `seednode::Int`.
	- Parameters (varied + constants): `b_cxn_peer::Float64`, `b_cxn_global::Float64`, `b_self::Float64`,
	`b_int::Float64`, `b_close::Float64`, `b_cls_x_smp::Float64`, `p_symp::Float64`, `beta::Float64`.
	- Outcomes per time step: `time::Float64`, `propeverinf::Float64`, `propcurrinf::Float64`, `proprec::Float64`.

	**Examples**
	```julia
	using Random, DataFrames

	Random.seed!(1234)
	df = replicate_sas_simulation(alst, vlst, 0.02, 14, 200, 20)

	# Aggregate to final timestep per (coef combo × itter)
	final_df = combine(groupby(df, [:coef_row, :itter])) do g
		g[argmax(g.time), :]
	end

	# Median curves over iterations at each design row (example: propeverinf)
	med_df = combine(groupby(df, [:coef_row, :time, :b_cxn_peer, :b_cxn_global, :b_self])) do g
		(propeverinf_med = median(g.propeverinf),
		propcurrinf_med = median(g.propcurrinf),
		proprec_med     = median(g.proprec))
	end
	first(med_df, 6)
	See Also
	compare_parameter_sweeps_probact_design_matrix, sirdif
	""" replicate_sas_simulation

#	Helper Function for simulation_comparer: read SAS CSV and normalize + aggregate
	function sas_simulation_aggregator(sas_csv_path::AbstractString)
		"""
		Args:
			sas_csv_path::AbstractString: File path to SAS output CSV with header:
				ITTER,BETA,TIME,PROPEVERINF,PROPCURRINF,PROPREC,case,
				B_INT,B_CXN_PEER,B_CXN_GLOBAL,B_CLOSE,B_SELF,B_CLS_X_SMP,SEEDNODE,maxinfected,ismax
		Returns:
			Tuple{DataFrame,DataFrame,DataFrame,DataFrame}:
				(df, med_time, final_per_run, med_final)

				- df::DataFrame
				Normalized long-form SAS table (all timesteps, all iterations).
				- med_time::DataFrame
				Median trajectories across iterations for each (time, b_cxn_peer, b_cxn_global, b_self).
				- final_per_run::DataFrame
				One final-timestep row per (b_cxn_peer, b_cxn_global, b_self, itter).
				- med_final::DataFrame
				Medians across iterations of the final-timestep outcomes per coefficient setting.
		Notes:
			- Keeps only columns needed for comparison/traceability; renames to Julia-style lowercase.
			- Assumes beta, b_int, b_close, b_cls_x_smp are constant within a coefficient group.
		"""

		#	Import
			df_raw = DataFrame(CSV.File(sas_csv_path))

		#	Select required columns (by SAS names)
			keep = [:ITTER, :BETA, :TIME, :PROPEVERINF, :PROPCURRINF, :PROPREC,
					:B_INT, :B_CXN_PEER, :B_CXN_GLOBAL, :B_CLOSE, :B_SELF, :B_CLS_X_SMP, :SEEDNODE]
			df = df_raw[!, keep]

		#	Rename to lowercase, Julia-style
			rename!(df, Dict(
				:ITTER=>:itter, :BETA=>:beta, :TIME=>:time,
				:PROPEVERINF=>:propeverinf, :PROPCURRINF=>:propcurrinf, :PROPREC=>:proprec,
				:B_INT=>:b_int, :B_CXN_PEER=>:b_cxn_peer, :B_CXN_GLOBAL=>:b_cxn_global,
				:B_CLOSE=>:b_close, :B_SELF=>:b_self, :B_CLS_X_SMP=>:b_cls_x_smp, :SEEDNODE=>:seednode
			))

		#	Enforce types (vectorized; no push!)
			df.itter        = Int.(df.itter)
			df.seednode     = Int.(df.seednode)

			df.beta         = Float64.(df.beta)
			df.time         = Float64.(df.time)
			df.propeverinf  = Float64.(df.propeverinf)
			df.propcurrinf  = Float64.(df.propcurrinf)
			df.proprec      = Float64.(df.proprec)
			df.b_int        = Float64.(df.b_int)
			df.b_cxn_peer   = Float64.(df.b_cxn_peer)
			df.b_cxn_global = Float64.(df.b_cxn_global)
			df.b_close      = Float64.(df.b_close)
			df.b_self       = Float64.(df.b_self)
			df.b_cls_x_smp  = Float64.(df.b_cls_x_smp)

		#	Median curves across iterations (per time & coef setting)
			g_time = groupby(df, [:time, :b_cxn_peer, :b_cxn_global, :b_self])
			med_time = combine(g_time,
				:beta        => first  => :beta,
				:b_int       => first  => :b_int,
				:b_close     => first  => :b_close,
				:b_cls_x_smp => first  => :b_cls_x_smp,
				:propeverinf => median => :propeverinf_med,
				:propcurrinf => median => :propcurrinf_med,
				:proprec     => median => :proprec_med
			)
			sort!(med_time, [:b_cxn_peer, :b_cxn_global, :b_self, :time])

		#	Final timestep per (coef combo × iteration)
			g_run = groupby(df, [:b_cxn_peer, :b_cxn_global, :b_self, :itter])
			final_per_run = combine(g_run) do gi
				i = argmax(gi.time)
				(time_final            = gi.time[i],
				propeverinf_final     = gi.propeverinf[i],
				propcurrinf_final     = gi.propcurrinf[i],
				proprec_final         = gi.proprec[i],
				beta                  = gi.beta[i],
				b_int                 = gi.b_int[i],
				b_close               = gi.b_close[i],
				b_cls_x_smp           = gi.b_cls_x_smp[i])
			end
			sort!(final_per_run, [:b_cxn_peer, :b_cxn_global, :b_self, :itter])

		#	Medians across iterations of those finals (per coef setting)
			g_coef = groupby(final_per_run, [:b_cxn_peer, :b_cxn_global, :b_self])
			med_final = combine(g_coef,
				:beta              => first  => :beta,
				:b_int             => first  => :b_int,
				:b_close           => first  => :b_close,
				:b_cls_x_smp       => first  => :b_cls_x_smp,
				:time_final        => median => :time_final_med,
				:propeverinf_final => median => :propeverinf_final_med,
				:propcurrinf_final => median => :propcurrinf_final_med,
				:proprec_final     => median => :proprec_final_med
			)
			sort!(med_final, [:b_cxn_peer, :b_cxn_global, :b_self])

		#	Return all requested tables
			return df, med_time, final_per_run, med_final
	end

#	Helper Function for simulation_comparer: Julia medians and finals
	function simulation_comparer_build_julia_medians(df::DataFrame)
		"""
		**Description**
		Aggregates the Julia simulation output (from `replicate_sas_simulation`) into:
		(1) median trajectories by time and coefficient setting,
		(2) final-timestep rows per run, and
		(3) medians of those finals per coefficient setting.

		**Usage**
		`simulation_comparer_build_julia_medians(df)`

		**Arguments**
		- `df::DataFrame`: Long-form table returned by `replicate_sas_simulation`.

		**Details**
		- `med_time`: groups by `(time, b_cxn_peer, b_cxn_global, b_self)` and takes medians of
		`propeverinf`, `propcurrinf`, `proprec` across iterations.
		- `final_per_run`: for each `(b_*, itter)` group, selects the row with the maximum `time`.
		- `med_final`: medians of the final-timestep outcomes across iterations per
		`(b_cxn_peer, b_cxn_global, b_self)`.

		**Value**
		Returns a `NamedTuple`:
		- `:med_time  => DataFrame` with columns
		`time, b_cxn_peer, b_cxn_global, b_self, propeverinf_med_julia, propcurrinf_med_julia, proprec_med_julia`.
		- `:final_per_run => DataFrame` with columns
		`b_cxn_peer, b_cxn_global, b_self, itter, time_final, propeverinf_final, propcurrinf_final, proprec_final`.
		- `:med_final => DataFrame` with columns
		`b_cxn_peer, b_cxn_global, b_self, time_final_med, propeverinf_final_med, propcurrinf_final_med, proprec_final_med`.

		**See Also**
		`sas_simulation_aggregator`, `replicate_sas_simulation`
		"""

		#	Median trajectories across iterations (per time & coef setting) ----
			g_time = groupby(df, [:time, :b_cxn_peer, :b_cxn_global, :b_self])
			med_time = combine(g_time,
				:propeverinf => median => :propeverinf_med_julia,
				:propcurrinf => median => :propcurrinf_med_julia,
				:proprec     => median => :proprec_med_julia
			)
			sort!(med_time, [:b_cxn_peer, :b_cxn_global, :b_self, :time])

		#	Final timestep per (coef combo × iteration) ----
		# 	Mirror SAS grouping keys for easy comparison
			g_run = groupby(df, [:b_cxn_peer, :b_cxn_global, :b_self, :itter])
			final_per_run = combine(g_run) do gi
				i = argmax(gi.time)
				(time_final        = gi.time[i],
				propeverinf_final = gi.propeverinf[i],
				propcurrinf_final = gi.propcurrinf[i],
				proprec_final     = gi.proprec[i])
			end
			sort!(final_per_run, [:b_cxn_peer, :b_cxn_global, :b_self, :itter])

		#	Medians of finals across iterations (per coef setting) ----
			g_coef = groupby(final_per_run, [:b_cxn_peer, :b_cxn_global, :b_self])
			med_final = combine(g_coef,
				:time_final        => median => :time_final_med,
				:propeverinf_final => median => :propeverinf_final_med,
				:propcurrinf_final => median => :propcurrinf_final_med,
				:proprec_final     => median => :proprec_final_med
			)
			sort!(med_final, [:b_cxn_peer, :b_cxn_global, :b_self])

		#	Return Median Datasets
			return (med_time = med_time, final_per_run = final_per_run, med_final = med_final)
	end

#	Main Wrapper: SAS ↔ Julia median-trajectory comparison + R plot
	function sas_simulation_comparer(
		alst, vlst,
		inf_r::Float64, rec_t::Int, maxtime::Int, num_iter::Int,
		sas_csv_path::AbstractString;
		outcome::Symbol = :propeverinf,
		display::String = ":100",
		save_path::Union{String,Nothing} = nothing,
		show_plot::Bool = true,
		seed::Union{Int,Nothing} = nothing)
		"""
		Args:
			alst, vlst: Network adjacency and weights
			inf_r: Base transmission probability
			rec_t: Recovery time in days
			maxtime: Maximum simulation days
			num_iter: Iterations per coefficient setting
			sas_csv_path: Path to SAS CSV output
		Keywords:
			outcome: :propeverinf | :propcurrinf | :proprec for plotting
			display: X11 display target
			save_path: PDF output path (nothing for X11)
			show_plot: Whether to open X11 window
			seed: RNG seed for reproducibility
		Returns:
			DataFrame with comparison table
		Notes:
			Changed to handle DataFrame infection_log from sirdif
		"""
		
		#	Seed (optional)
			if seed !== nothing
				Random.seed!(seed)
			end

		#	Read SAS CSV + aggregates
			df_sas, med_time_sas, _final_per_run_sas, _med_final_sas = sas_simulation_aggregator(sas_csv_path)

		#	Matching parameters from SAS file
			num_iter_sas  = length(unique(df_sas.itter))
			maxtime_sas   = maximum(df_sas.time)

		#	Run Julia simulation (now returns DataFrame)
			df_julia = replicate_sas_simulation(alst, vlst, inf_r, rec_t, Int(maxtime_sas), num_iter_sas)

		#	Aggregate Julia medians/finals
			julia_summ     = simulation_comparer_build_julia_medians(df_julia)
			med_time_julia = julia_summ.med_time

		#	Harmonize SAS median column names → *_med_sas
			rename!(med_time_sas, Dict(
				:propeverinf_med => :propeverinf_med_sas,
				:propcurrinf_med => :propcurrinf_med_sas,
				:proprec_med     => :proprec_med_sas
			))

		#	Validate coefficients only (not time grid)
			coef_keys    = [:b_cxn_peer, :b_cxn_global, :b_self]
			levels_sas   = sort(unique(eachrow(med_time_sas[:, coef_keys])))
			levels_julia = sort(unique(eachrow(med_time_julia[:, coef_keys])))
			@assert levels_sas == levels_julia "Coefficient sets differ."

		#	Per-time counts contributing to medians (SAS & Julia)
			cnt_sas = combine(groupby(df_sas, [:time, :b_cxn_peer, :b_cxn_global, :b_self]),
				nrow => :n_sas)
			cnt_jul = combine(groupby(df_julia, [:time, :b_cxn_peer, :b_cxn_global, :b_self]),
				nrow => :n_julia)

		#	Join on (time + coefficients) using intersection of times
			comp = innerjoin(
				med_time_sas, med_time_julia,
				on = [:time, :b_cxn_peer, :b_cxn_global, :b_self],
				makeunique = true
			)

		#	Attach counts to the comparison rows
			comp = innerjoin(comp, cnt_sas, on = [:time, :b_cxn_peer, :b_cxn_global, :b_self])
			comp = innerjoin(comp, cnt_jul, on = [:time, :b_cxn_peer, :b_cxn_global, :b_self])

		#	Constants used in Julia run (held fixed)
			const_p_symp      = 0.75
			const_b_int       = -0.1
			const_b_close     =  1.0
			const_b_cls_x_smp = -0.1

		#	Absolute deltas: abs(SAS − Julia)
			propeverinf_delta = abs.(comp.propeverinf_med_sas .- comp.propeverinf_med_julia)
			propcurrinf_delta = abs.(comp.propcurrinf_med_sas .- comp.propcurrinf_med_julia)
			proprec_delta     = abs.(comp.proprec_med_sas     .- comp.proprec_med_julia)

		#	Assemble comparison table
			nrows     = nrow(comp)
			beta_v    = fill(inf_r, nrows)
			p_symp_v  = fill(const_p_symp, nrows)
			b_int_v   = fill(const_b_int, nrows)
			b_close_v = fill(const_b_close, nrows)
			b_cls_v   = fill(const_b_cls_x_smp, nrows)
			rec_t_v   = fill(rec_t, nrows)
			maxtime_v = fill(maxtime, nrows)

			comp_table = DataFrame(
				beta        = beta_v,
				p_symp      = p_symp_v,
				b_int       = b_int_v,
				b_close     = b_close_v,
				b_cls_x_smp = b_cls_v,
				rec_t       = rec_t_v,
				maxtime     = maxtime_v,

				b_cxn_peer   = comp.b_cxn_peer,
				b_cxn_global = comp.b_cxn_global,
				b_self       = comp.b_self,
				time         = comp.time,

				propeverinf_med_sas = comp.propeverinf_med_sas,
				propcurrinf_med_sas = comp.propcurrinf_med_sas,
				proprec_med_sas     = comp.proprec_med_sas,

				propeverinf_med_julia = comp.propeverinf_med_julia,
				propcurrinf_med_julia = comp.propcurrinf_med_julia,
				proprec_med_julia     = comp.proprec_med_julia,

				propeverinf_delta = propeverinf_delta,
				propcurrinf_delta = propcurrinf_delta,
				proprec_delta     = proprec_delta,

				n_sas   = comp.n_sas,
				n_julia = comp.n_julia
			)

		#	Prep labels for plotting
			@rput comp_table display
			ycol_j = outcome === :propeverinf ? "propeverinf_med_julia" :
					outcome === :propcurrinf ? "propcurrinf_med_julia" :
											"proprec_med_julia"
			ycol_s = outcome === :propeverinf ? "propeverinf_med_sas" :
					outcome === :propcurrinf ? "propcurrinf_med_sas" :
											"proprec_med_sas"
			y_label = outcome === :propeverinf ? "Proportion Ever Infected" :
					outcome === :propcurrinf ? "Proportion Currently Infected" :
												"Proportion Recovered"

			title_label    = "SAS vs Julia Median Trajectories"
			subtitle_label = "Columns: Peer effect   Rows: Global effect   Lines: Symptomatic"
			param_text = "β=$(round(inf_r, digits=4))  p_symp=$(const_p_symp)  b_int=$(const_b_int)  " *
						"b_close=$(const_b_close)  b_cls_x_smp=$(const_b_cls_x_smp)  rec_t=$(rec_t)  maxtime=$(maxtime)"

			@rput comp_table ycol_j ycol_s y_label title_label subtitle_label param_text save_path

		#	Open graphics device if saving / showing
			if save_path !== nothing
				R"""
				pdf(save_path, width=10, height=12)
				"""
			elseif show_plot
				R"""
				Sys.setenv(DISPLAY=':100')
				Sys.getenv("DISPLAY")
				X11(type = "cairo", width=10, height=12)
				"""
			end

		#	R plotting (solid = Julia, dashed = SAS)
			R"""
			col_j <- ycol_j
			col_s <- ycol_s
			ylab  <- y_label
			ttl   <- title_label
			subt  <- subtitle_label
			param <- param_text
			spath <- save_path

			df <- comp_table
			df[["b_cxn_peer"]]   <- as.numeric(df[["b_cxn_peer"]])
			df[["b_cxn_global"]] <- as.numeric(df[["b_cxn_global"]])
			df[["b_self"]]       <- as.numeric(df[["b_self"]])
			df[["time"]]         <- as.numeric(df[["time"]])

			peer_levels <- sort(unique(df[["b_cxn_peer"]]))
			glob_levels <- sort(unique(df[["b_cxn_global"]]))
			self_levels <- sort(unique(df[["b_self"]]))

			layout(matrix(1:(length(glob_levels)*length(peer_levels)),
						nrow=length(glob_levels), byrow=TRUE))
			par(family="serif", mar=c(3.5,4,2.5,1.5), las=1)

			dataplot_tick_function <- function(major_tick_length=0.035, minor_tick_ratio=0.25){
				packages <- c('Hmisc')
				missing <- setdiff(packages, rownames(installed.packages()))
				if (length(missing) > 0) install.packages(missing, repos = "https://cloud.r-project.org")
				Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = minor_tick_ratio)
				Hmisc::minor.tick(nx = 2, ny = 2, tick.ratio = -minor_tick_ratio)
				axis(2, tck=1, tck=-major_tick_length, labels = FALSE)
				axis(1, tck=1, tck=-major_tick_length, labels = FALSE)
			}

			pal_base <- c("#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#4477AA", "#BBBBBB")
			if (length(self_levels) <= length(pal_base)) {
				pal <- pal_base[seq_along(self_levels)]
			} else {
				pal <- grDevices::colorRampPalette(pal_base)(length(self_levels))
			}

			for (gi in seq_along(glob_levels)) {
				gval <- glob_levels[gi]
				for (pi in seq_along(peer_levels)) {
					pval <- peer_levels[pi]
					sub  <- df[df[["b_cxn_global"]] == gval & df[["b_cxn_peer"]] == pval, ]

					ylim <- range(c(sub[[ col_j ]], sub[[ col_s ]]), na.rm=TRUE)
					ylim[1] <- min(0, ylim[1])

					plot(NA, type="n",
						xlim=range(sub[["time"]], na.rm=TRUE), ylim=ylim,
						xlab=" ", ylab=ylab, tck=0.015, xaxt='n', bty='L', las=1,
						main=paste0("Peer effect = ", sprintf("%.2f", pval),
									"   |   Global effect = ", sprintf("%.2f", gval)))
					mtext(side = 1, text = 'Time', col = "black", line = 2.45, cex = 0.75, family='serif')
					axis(1, padj=0.75, tck=0.015)
					dataplot_tick_function()

					for (si in seq_along(self_levels)) {
						sval <- self_levels[si]
						ss <- sub[sub[["b_self"]] == sval, ]
						ss <- ss[order(ss[["time"]]), ]

						lines(ss[["time"]], ss[[ col_j ]], col=pal[si], lwd=2, lty=1)
						lines(ss[["time"]], ss[[ col_s ]], col=pal[si], lwd=2, lty=2)
					}
				}
			}

			if (!is.null(spath)) dev.off()
			"""

		#	Return comparison table
			return comp_table
	end
	@doc raw"""
	**Description**
	Compares SAS and Julia simulations on a fixed network by running the Julia
	design-driven simulation, reading the SAS CSV, and joining **per-time median
	trajectories** on shared coefficient keys. Produces an R panel plot (X11/PDF)
	and returns a single table with constants, SAS medians, Julia medians, and
	absolute deltas.

	**Usage**
	`sas_simulation_comparer(alst, vlst, inf_r, rec_t, maxtime, num_iter, sas_csv_path;
							outcome=:propeverinf, display=":100",
							save_path=nothing, show_plot=true, seed=nothing)`

	**Arguments**
	- `alst`, `vlst`: Network adjacency and weights (matrix or vector-of-vectors; same form).
	- `inf_r::Float64`: Base transmission probability per contact (β).
	- `rec_t::Int`: Recovery time in days.
	- `maxtime::Int`: Maximum simulation days (≥ 1).
	- `num_iter::Int`: Iterations per coefficient setting (≥ 1).
	- `sas_csv_path::AbstractString`: Path to SAS PROC IML output CSV.

	**Keywords**
	- `outcome::Symbol`: `:propeverinf` (default), `:propcurrinf`, or `:proprec` for plotting.
	- `display::String`: X11 display target for R (e.g., `":100"`).
	- `save_path::Union{String,Nothing}`: PDF output path; if set, no X11 window opens.
	- `show_plot::Bool`: Whether to open an X11 window when `save_path` is `nothing`.
	- `seed::Union{Int,Nothing}`: If set, seeds Julia RNG for reproducibility.

	**Details**
	- Validates exact match of `(time, b_cxn_peer, b_cxn_global, b_self)` grids.
	- Faceting: rows = **Global effect** (`b_cxn_global`), cols = **Peer effect** (`b_cxn_peer`).
	- Line color encodes **Symptomatic** (`b_self`) level; Julia = solid, SAS = dashed.
	- Bottom annotation shows constants: `β, p_symp, b_int, b_close, b_cls_x_smp, rec_t, maxtime`.

	**Value**
	Returns a `DataFrame` with columns (left→right: constants → coefficients → SAS → Julia → deltas):
	- `beta, p_symp, b_int, b_close, b_cls_x_smp, rec_t, maxtime`
	- `b_cxn_peer, b_cxn_global, b_self, time`
	- `propeverinf_med_sas, propcurrinf_med_sas, proprec_med_sas`
	- `propeverinf_med_julia, propcurrinf_med_julia, proprec_med_julia`
	- `propeverinf_delta, propcurrinf_delta, proprec_delta` (absolute deltas).

	**See Also**
	`replicate_sas_simulation`, `simulation_comparer_build_julia_medians`, `sas_simulation_aggregator`
	""" sas_simulation_comparer

######################
#   FUNCTION TESTS   #
######################

#   arithmetic_mode 

#   Test 1: Simple Vector (PASSED)
    test1 = [1, 2, 3, 2, 2, 4, 5]
    result1 = arithmetic_mode(test1)
    println("Test 1 - Input: $test1")
    println("Mode: $result1")
    println("Expected: 2")
    println()

#   Test 2: Array with missing values (keeping missing)
    test2 = [1, 2, missing, 2, 3, missing, missing]
    result2 = arithmetic_mode(test2, na_rm=false)
    println("Test 2 - Input with missing (na_rm=false): $test2")
    println("Mode: $result2")
    println("Expected: missing (most frequent)")
    println()

#   Test 3: Array with missing values (removing missing)
    test3 = [1, 2, missing, 2, 3, missing, 2]
    result3 = arithmetic_mode(test3, na_rm=true)
    println("Test 3 - Input with missing (na_rm=true): $test3")
    println("Mode: $result3")
    println("Expected: 2")
    println()

#   Test 4: Float array
    test4 = [1.5, 2.3, 1.5, 3.7, 1.5, 2.3]
    result4 = arithmetic_mode(test4)
    println("Test 4 - Float input: $test4")
    println("Mode: $result4")
    println("Expected: 1.5")
    println()

#   Test 5: String array
    test5 = ["apple", "banana", "apple", "cherry", "apple", "banana"]
    result5 = arithmetic_mode(test5)
    println("Test 5 - String input: $test5")
    println("Mode: $result5")
    println("Expected: apple")

#   color_assignment

#   Test 1: Simple test with values across the range
    test1 = [0.02, 0.07, 0.12, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43, 0.48, 0.52, 0.58]
    result1 = color_assignment(test1)
 

    println("Test 1 - Input values: $test1")
    println("Colors (optimized): $result1")
    println()

#   Test 2: Edge cases
    test2 = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
    result2 = color_assignment(test2)
    println("Test 2 - Edge values: $test2")
    println("Colors: $result2")
    println()

#   Test 3: Random values
    test3 = rand(10) * 0.6  # Random values between 0 and 0.6
    result3 = color_assignment(test3)
    println("Test 3 - Random values: $(round.(test3, digits=3))")
    println("Colors: $result3")
    println()

#   Test 4: Single value
    test4 = [0.27]
    result4 = color_assignment(test4)
    println("Test 4 - Single value: $test4")
    println("Color (optimized): $result4")
    println("Expected: navajowhite (for range [0.25, 0.30))")
    println()

#   Parameter Space Exploration

#   Full test with output
    test_results = test_sim_parm_space(verbose=true, export_csv=false)

#   Plot results
    results = sim_parm_space(-0.1, -1.5, -3.5, -0.5, 1.0, -0.1, 3, 5, 500)
    plot_parameter_space(results, edgemax=3)

#   Simulation Preparation: sim_prep()

#   Load network from GraphML file
	network_data = sim_prep("/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml")
	alst = network_data.alst
	vlst = network_data.vlst
	
#   Check network structure
	n_nodes = size(alst, 1)
	max_degree = size(alst, 2) - 1
	println("Network has $n_nodes nodes with max degree $max_degree")

#   Simulation Execution: sirdif()

#   SIR Simulation Tests
    sir_results = test_sir_model(network_data; seed=123, verbose=true)
   
	let log = sir_results["results"]["infection_log"]
		#   Daily growth factors up to peak
            peak_idx = findfirst(==(maximum(log[:,2])), log[:,2])
            growth = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
            println("Median pre-peak growth factor: ", round(median(growth), digits=3))

		#   Effective attack rate by peak day
            prop_ever_at_peak = log[peak_idx, 3]
            println("Prop. ever infected by peak day: ", round(prop_ever_at_peak, digits=3))
	end
    
    let log = sir_results["results"]["infection_log"]; ni = log[:,2]
        println("NI (first 30 days): ", join(round.(ni[1:min(end,30)]), ", "))
    end

    let log = sir_results["results"]["infection_log"]
        peak_idx = findfirst(==(maximum(log[:,2])), log[:,2])
        growth   = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
        println("Median pre-peak growth: ", round(median(growth), digits=3))
        println("Prop. ever infected at peak: ", round(log[peak_idx, 3], digits=3))
        println("Peak day: ", peak_idx, "   Peak NI: ", log[peak_idx, 2])
    end

#	Small Speed Test
	directory = "/workspace/data/sim_test_data/synthetic_village_m_graphml"
	graphml_files = readdir(directory)
	results_list = Vector{Dict}(undef, length(graphml_files))
	iteration_number = 0
	@showprogress 1 "Running simulations…" for i in eachindex(graphml_files)
		# 	Importing Network
			iteration_number = i
			network_data = sim_prep(string(directory, "/", graphml_files[i]))

		#	Use matrix format directly from sim_prep (sirdif will convert internally)
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)

		#	Set simulation parameters
			params = Dict(
				"inf_r"        => 0.33,
				"rec_t"        => 14,
				"maxtime"      => 200,
				"p_symp"       => 0.5,
				"b_int"        => -0.1,
				"b_close"      => 1.0,
				"b_cxn_peers"  => -0.5,
				"b_cxn_total"  => -3.5,
				"b_cxn_symp"   => -1.5,
				"b_cls_x_smp"  => -0.1
			)

		#	Fixed seed nodes (within bounds)
			raw_seeds = [5, 10, 15, 20, 25]
			seed_nodes = [s for s in raw_seeds if 1 ≤ s ≤ n]
			isempty(seed_nodes) && error("Provided seed nodes $(raw_seeds) are out of bounds for n = $n.")

		#	Run simulation (SIR: permanent immunity)
			results = sirdif(
				alst, vlst, seed_nodes,
				params["inf_r"], params["rec_t"], params["maxtime"], params["p_symp"],
				params["b_int"], params["b_close"], params["b_cxn_peers"],
				params["b_cxn_total"], params["b_cxn_symp"], params["b_cls_x_smp"];
				transmission_method = :weighted,
				immunity_duration   = nothing
			)

		#	Extract results
			infection_log = results["infection_log"]
			runtime = results["total_time"]

		#	Key metrics
			final_row     = infection_log[end, :]
			final_size    = final_row[3]
			peak_infected = maximum(infection_log[:, 2])
			peak_time     = findfirst(==(peak_infected), infection_log[:, 2])

		#	Populate Results Vector
			results_list[i] = Dict(
				"results" => results,
				"metrics" => Dict(
					"final_size"    => final_size,
					"peak_infected" => peak_infected,
					"peak_time"     => peak_time,
					"duration"      => size(infection_log, 1) - 1
				),
				"params" => params,
				"seeds"  => seed_nodes
			)
	end	

#	Examing Outputs
	network_26_results = results_list[1000]

	let log = network_26_results["results"]["infection_log"]
		#   Daily growth factors up to peak
            peak_idx = findfirst(==(maximum(log[:,2])), log[:,2])
            growth = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
            println("Median pre-peak growth factor: ", round(median(growth), digits=3))

		#   Effective attack rate by peak day
            prop_ever_at_peak = log[peak_idx, 3]
            println("Prop. ever infected by peak day: ", round(prop_ever_at_peak, digits=3))
	end

    
    let log = network_26_results["results"]["infection_log"]; ni = log[:,2]
        println("NI (first 30 days): ", join(round.(ni[1:min(end,30)]), ", "))
    end

    let log = network_26_results["results"]["infection_log"]
        peak_idx = findfirst(==(maximum(log[:,2])), log[:,2])
        growth   = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
        println("Median pre-peak growth: ", round(median(growth), digits=3))
        println("Prop. ever infected at peak: ", round(log[peak_idx, 3], digits=3))
        println("Peak day: ", peak_idx, "   Peak NI: ", log[peak_idx, 2])
    end

#   SIRS Simulation Tests
#   sirs_results = test_sirs_model(network_data; immunity_days=14, seed=42, verbose=true)

#	Pre-peak growth + attack rate at peak
	let log = sirs_results["results"]["infection_log"]
		#	Locate peak (by NI)
			peak_idx = findfirst(==(maximum(log[:, 2])), log[:, 2])

		#	Daily growth factors up to peak
			growth = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
			med_growth = round(median(growth), digits=3)

		#	Effective attack rate by peak day (cumulative)
			prop_ever_at_peak = round(log[peak_idx, 3], digits=3)

		println("Median pre-peak growth factor: ", med_growth)
		println("Prop. ever infected by peak day: ", prop_ever_at_peak)
	end

#	Quick glance at early incidence (NI)
	let log = sirs_results["results"]["infection_log"]; ni = log[:, 2]
		println("NI (first 30 days): ", join(round.(ni[1:min(end, 30)]), ", "))
	end

#	Wave detection (smoothed), timing vs. recovery+immunity gap
	let log = sirs_results["results"]["infection_log"]
        #   Specifying number infected
		    ni = log[:, 2]

		#	Smooth lightly (5-day moving average) to suppress single-day noise
			smoothed = similar(ni)
			window = 5
			halfwin = fld(window, 2)
			for i in eachindex(ni)
				l = max(1, i - halfwin)
				r = min(length(ni), i + halfwin)
				smoothed[i] = sum(ni[l:r]) / (r - l + 1)
			end

		#	Peak selection rules (match test logic)
			immunity_days = 30					# set to your call’s value
			rec_t = 14							# from R settings
			min_gap = immunity_days + rec_t
			peak_height_thresh = max(5.0, 0.05 * maximum(smoothed))

		#	Candidate local maxima on smoothed series
			candidates = Int[]
			for i in 2:(length(smoothed)-1)
				if (smoothed[i] > smoothed[i-1]) && (smoothed[i] > smoothed[i+1]) && (smoothed[i] ≥ peak_height_thresh)
					push!(candidates, i)
				end
			end

		#	Enforce minimum gap: keep taller peak if two are too close
			peaks = Int[]
			for idx in candidates
				if isempty(peaks) || (idx - peaks[end] ≥ min_gap)
					push!(peaks, idx)
				else
					if smoothed[idx] > smoothed[peaks[end]]
						peaks[end] = idx
					end
				end
			end

		#	Print wave summary
			if isempty(peaks)
				println("Waves: none detected")
			else
				heights = round.(smoothed[peaks]; digits=1)
				println("Waves detected (day → smoothed NI):")
				for (d, h) in zip(peaks, heights)
					println("  day ", d, " → ", h)
				end
				if length(peaks) > 1
					gaps = diff(peaks)
					println("Inter-wave gaps (days): ", join(gaps, ", "))
					println("All gaps ≥ ", min_gap, "? ", all(gaps .≥ min_gap) ? "YES" : "NO")
				end
			end
	end

#	Reinfection signal (from cumulative ever-infected measure)
#	Note: in SIRS, prop_ever_inf = (NI + NR)/n can exceed 1.0 if people are infected multiple times.
	let log = sirs_results["results"]["infection_log"]
        #   Specifying Parameters
            n = size(sirs_results["results"]["infection_log"], 1)	# rows in log, not population size
            final_prop_ever = log[end, 3]
            has_reinfection = final_prop_ever > 1.0
            println("Reinfection indicated (prop_ever_inf > 1 at end): ", has_reinfection ? "YES" : "NO")

		#	Estimate total infections and reinfections (uses population from params already printed earlier)
		#	You can pass population size directly if you prefer: pop_n = size(network_data.alst, 1)
            pop_n = size(network_data.alst, 1)
            total_infections_est = round(Int, final_prop_ever * pop_n)
            est_reinfections = max(0, total_infections_est - pop_n)
            println("Total infections (est.): ", total_infections_est, "   Reinfections (est.): ", est_reinfections)

		#	First day cumulative ever-infected exceeds population (if any)
            exceed_idx = findfirst(x -> x > 1.0, log[:, 3])
            if exceed_idx === nothing
                println("Day crossing >100% cumulative: none within horizon")
            else
                println("Day crossing >100% cumulative: day ", exceed_idx)
            end
	end

#	Peak recap (SIRS)
	let log = sirs_results["results"]["infection_log"]
		peak_idx = findfirst(==(maximum(log[:, 2])), log[:, 2])
		growth   = log[2:peak_idx, 2] ./ max.(1, log[1:peak_idx-1, 2])
		println("Median pre-peak growth: ", round(median(growth), digits=3))
		println("Prop. ever infected at peak: ", round(log[peak_idx, 3], digits=3))
		println("Peak day: ", peak_idx, "   Peak NI: ", log[peak_idx, 2])
	end

#   Comparisative Tests: Replication Analysis

#   Load network from GraphML file
	network_data = sim_prep("/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml", sas_transformation=true)
	alst = network_data.alst
	vlst = network_data.vlst

#	Parameters (SAS IML-like) ---
	Random.seed!(1234)
	inf_r   = 0.02      # transmission prob (bt)
	rec_t   = 14        # days to recovery
	maxtime = 200       # simulation horizon
	num_iter = 200       # <- sanity run
	sas_csv = "/workspace/data/sim_test_data/SAS_Simulation_Outputs.csv"

#	Executing Comparison Simulations
	comp_table = sas_simulation_comparer(
		alst, vlst,                      # your fixed network
		inf_r, rec_t, maxtime, num_iter, # replication knobs
		sas_csv;                         # SAS results
		outcome = :propeverinf,          # or :propcurrinf / :proprec
		seed    = 1234,
		show_plot = true,
		display   = ":100",              # X11 target
		save_path = nothing
	)

#	Looking at Maximum Difference
	maximum(comp_table.propeverinf_delta)

########################
#   DIAGNOSTIC TESTS   #
########################

#	Rare Events Analysis
	function rare_event_diagnostic(
			comp::DataFrame;
			outcome::Symbol = :propeverinf,		# :propeverinf | :propcurrinf | :proprec
			count_metric::Symbol = :n_sas,		# :n_sas | :n_julia | :n_min | :n_hmean
			bins::Int = 5)
		#	Pick the delta column implied by `outcome`
			delta_col = Symbol(string(outcome), "_delta")
			@assert hasproperty(comp, delta_col) "Column $(delta_col) not found in comp table."

		#	Ensure count columns are present
			@assert hasproperty(comp, :n_sas)   "comp missing :n_sas"
			@assert hasproperty(comp, :n_julia) "comp missing :n_julia"

		#	Select the count vector per requested metric
			ns = comp.n_sas
			nj = comp.n_julia

			counts = if count_metric === :n_sas
				ns
			elseif count_metric === :n_julia
				nj
			elseif count_metric === :n_min
				min.(ns, nj)
			elseif count_metric === :n_hmean
				den = ns .+ nj
				h   = similar(den, Float64)
				@inbounds for i in eachindex(den)
					h[i] = (den[i] == 0) ? 0.0 : (2.0 * ns[i] * nj[i]) / den[i]
				end
				h
			else
				throw(ArgumentError("Unsupported count_metric=$(count_metric)"))
			end

		#	Filter to valid rows (finite numbers, nonnegative counts)
			delta = comp[!, delta_col]
			ok = trues(length(delta))
			@inbounds for i in eachindex(ok)
				ok[i] &= isfinite(delta[i]) & isfinite(counts[i]) & (counts[i] ≥ 0)
			end
			if !any(ok)
				return (note = "No valid rows to test.", spearman = nothing, by_bin = DataFrame())
			end

		#	Slice, convert to Float64
			d = Float64.(delta[ok])
			c = Float64.(counts[ok])

		#	Spearman correlations (StatsBase.corspearman)
			ρ_count = StatsBase.corspearman(d, c)

			invc = similar(c)
			@inbounds for i in eachindex(invc)
				invc[i] = 1.0 / max(c[i], 1.0)		# treat 0 as 1 to avoid Inf
			end
			ρ_inv = StatsBase.corspearman(d, invc)

		#	Scarcity bins (equal-frequency “ntiles” by counts)
			n   = length(c)
			idx = sortperm(c, lt = <)
			tile = zeros(Int, n)
			@inbounds for b in 1:bins
				i_lo = floor(Int, (b - 1) * n / bins) + 1
				i_hi = floor(Int, b * n / bins)
				for j in i_lo:i_hi
					tile[idx[j]] = b
				end
			end

			tdf = DataFrame(bin = tile, counts = c, delta = d)
			by_bin = combine(groupby(tdf, :bin),
				:counts => minimum  => :count_min,
				:counts => median   => :count_med,
				:counts => maximum  => :count_max,
				:delta  => length   => :n_rows,
				:delta  => median   => :delta_median,
				:delta  => mean     => :delta_mean,
				:delta  => (x -> quantile(x, 0.90)) => :delta_p90,
				:delta  => maximum  => :delta_max,
			)
			sort!(by_bin, :bin)

		#	Return concise summary
			return (
				outcome      = outcome,
				count_metric = count_metric,
				rows_used    = sum(by_bin.n_rows),
				spearman     = (
					rho_delta_vs_count     = ρ_count,	# expect negative if “more data ⇒ smaller delta”
					rho_delta_vs_inv_count = ρ_inv		# expect positive if “rarer ⇒ bigger delta”
				),
				by_bin = by_bin
			)
	end
	
	diag = rare_event_diagnostic(comp_table; outcome=:propeverinf, count_metric=:n_sas, bins=5)
	println(diag.spearman)
	first(diag.by_bin, 5)

#	Within-panel rarity diagnostics (Spearman)
	function rarity_within_panel_diagnostics(
		comp_table::DataFrame;
		outcome::Symbol = :propeverinf,
		count_strategy::Symbol = :min,
		return_overall::Bool = true)
		
		#	Validate input
			required_cols = [:b_cxn_peer, :b_cxn_global, :n_sas, :n_julia]
			missing_cols = setdiff(required_cols, Symbol.(names(comp_table)))
			if !isempty(missing_cols)
				error("Missing required columns: $(missing_cols)")
			end

		#	Select delta column based on outcome
			delta_col = Symbol("$(outcome)_delta")
			if !(String(delta_col) in names(comp_table))
				error("Delta column $(delta_col) not found")
			end

		#	Extract data
			delta_all = Float64.(comp_table[!, delta_col])
			n_sas = Int.(comp_table.n_sas)
			n_julia = Int.(comp_table.n_julia)

		#	Select count series
			count_series = if count_strategy === :min
				min.(n_sas, n_julia)
			elseif count_strategy === :sas
				n_sas
			elseif count_strategy === :julia
				n_julia
			else
				error("Invalid count_strategy: $(count_strategy)")
			end

		#	Filter valid rows
			valid_mask = .!ismissing.(delta_all) .& .!ismissing.(count_series)
			delta_use = delta_all[valid_mask]
			count_use = Float64.(count_series[valid_mask])

		#	Overall correlation
			overall = if return_overall && length(delta_use) >= 2
				(
					rho_delta_vs_count = corspearman(delta_use, count_use),
					rho_delta_vs_invcount = corspearman(delta_use, 1.0 ./ max.(count_use, 1.0)),
					n = length(delta_use)
				)
			else
				(rho_delta_vs_count = NaN, rho_delta_vs_invcount = NaN, n = length(delta_use))
			end

		#	Per-panel analysis
			grouped = groupby(comp_table, [:b_cxn_peer, :b_cxn_global])
			per_panel = combine(grouped) do gi
				#	Extract panel data
					gi_delta = Float64.(gi[!, delta_col])
					gi_ns = Int.(gi.n_sas)
					gi_nj = Int.(gi.n_julia)
					
				#	Panel counts
					gi_count = if count_strategy === :min
						min.(gi_ns, gi_nj)
					elseif count_strategy === :sas
						gi_ns
					else
						gi_nj
					end

				#	Valid data
					valid = .!ismissing.(gi_delta) .& .!ismissing.(gi_count)
					if !any(valid)
						return DataFrame(
							rho_delta_vs_count = NaN,
							rho_delta_vs_invcount = NaN,
							n = 0
						)
					end

					dp = gi_delta[valid]
					cp = Float64.(gi_count[valid])

				#	Calculate correlations
					if length(dp) >= 2 && length(unique(dp)) >= 2 && length(unique(cp)) >= 2
						DataFrame(
							rho_delta_vs_count = corspearman(dp, cp),
							rho_delta_vs_invcount = corspearman(dp, 1.0 ./ max.(cp, 1.0)),
							n = sum(valid)
						)
					else
						DataFrame(
							rho_delta_vs_count = NaN,
							rho_delta_vs_invcount = NaN,
							n = sum(valid)
						)
					end
			end

		#	Return results
			return (overall = overall, per_panel = per_panel)
	end

#	Summarize the rarity effect patterns
	function summarize_rarity_effects(per_panel::DataFrame)
		#	Add effect strength categories
			per_panel.peer_strength = abs.(per_panel.b_cxn_peer)
			per_panel.global_strength = abs.(per_panel.b_cxn_global)
			per_panel.combined_strength = per_panel.peer_strength .+ per_panel.global_strength
			
		#	Key findings
			println("RARITY EFFECT SUMMARY")
			println("=" ^ 60)
			
		#	Overall pattern
			mean_rho = mean(per_panel.rho_delta_vs_invcount)
			println("Overall mean ρ(delta, 1/count): ", round(mean_rho, digits=3))
			println("  → Positive correlation confirms: Lower counts → Larger errors\n")
			
		#	Most affected panels
			println("Most Affected Panels (strongest rarity effect):")
			sorted = sort(per_panel, :rho_delta_vs_invcount, rev=true)
			for i in 1:3
				row = sorted[i, :]
				println("  $(i). Peer=$(row.b_cxn_peer), Global=$(row.b_cxn_global)")
				println("     ρ = $(round(row.rho_delta_vs_invcount, digits=3))")
			end
			
		#	Least affected panels  
			println("\nLeast Affected Panels:")
			for i in 1:3
				row = sorted[end-i+1, :]
				println("  $(i). Peer=$(row.b_cxn_peer), Global=$(row.b_cxn_global)")
				println("     ρ = $(round(row.rho_delta_vs_invcount, digits=3))")
			end
			
		#	Pattern analysis
			println("\nPattern Analysis:")
		
		#	By peer effect
			peer_groups = combine(groupby(per_panel, :b_cxn_peer), 
				:rho_delta_vs_invcount => mean => :mean_rho)
			sort!(peer_groups, :b_cxn_peer)
			println("  By Peer Effect:")
			for row in eachrow(peer_groups)
				println("    Peer=$(row.b_cxn_peer): mean ρ = $(round(row.mean_rho, digits=3))")
			end
			
		#	Correlation with effect strengths
			println("\nEffect Strength Correlations:")
			println("  ρ(rarity effect, |peer effect|): ", 
				round(cor(per_panel.rho_delta_vs_invcount, per_panel.peer_strength), digits=3))
			println("  ρ(rarity effect, |global effect|): ", 
				round(cor(per_panel.rho_delta_vs_invcount, per_panel.global_strength), digits=3))
			println("  ρ(rarity effect, combined |effects|): ", 
				round(cor(per_panel.rho_delta_vs_invcount, per_panel.combined_strength), digits=3))
			
			return per_panel
	end

# 	comp_table is what sas_simulation_comparer returns
	res = rarity_within_panel_diagnostics(comp_table; outcome=:propeverinf, count_strategy=:min)
	res.overall

	enhanced_panel = summarize_rarity_effects(res.per_panel)

#	Trace Analysis

#	Load SAS Trace CSV
	function load_sas_trace(sas_log::AbstractString)::DataFrame
		"""
		Args:
			sas_log::AbstractString: Path to SIR_SocialFeedback_Infection_Log.csv
		Returns:
			DataFrame: All attempt-level rows with correct column types
		Notes:
			Requires CSV.jl and DataFrames.jl.
		"""
		#	Read the CSV.
			df = CSV.read(sas_log, DataFrame; normalizenames=true)

		#	Coerce expected numeric columns (robust to SAS formatting).
			numcols = [:IT,:T,:J,:EGO,:K,:COLPOS,:ALTER,:EDGWGT,:PEERSINF,:NINF,
			           :ISSYMPTOMATIC,:LWACT,:PROB_ACT,:ACTIVE,:TRANPROB,:TRANS,
			           :BECAME_INFECTED,:NBRS]
			for c in numcols
				if c ∈ names(df)
					df[!, c] = parse.(Float64, string.(df[!, c]))
				end
			end

		#	Cast integer-like columns to Int (no rounding, just trunc).
			intcols = [:IT,:T,:J,:EGO,:K,:COLPOS,:ALTER,:PEERSINF,:NINF,
			           :ISSYMPTOMATIC,:ACTIVE,:TRANS,:BECAME_INFECTED,:NBRS]
			for c in intcols
				if c ∈ names(df)
					df[!, c] = Int.(df[!, c])
				end
			end

		#	Return.
			return df
	end

	sas_log = "/workspace/data/sim_test_data/SIR_SocialFeedback_Infection_Log.csv"
	df = load_sas_trace(sas_log)
	size(df), names(df)

#	Basic Sanity Checks
	function trace_sanity_checks(df; beta::Float64, atol::Float64=1e-9)
		# 	Expected PROB_ACT from LWACT
			expected_prob = 1.0 ./ (1.0 .+ exp.(-df.LWACT))
			prob_err = abs.(df.PROB_ACT .- expected_prob)
			max_prob_err = maximum(skipmissing(prob_err))
			ok_prob = max_prob_err ≤ atol

		# 	Expected TRANPROB from beta and EDGWGT
			expected_tranprob = 1 .- (1 .- beta) .^ df.EDGWGT
			tran_err = abs.(df.TRANPROB .- expected_tranprob)
			max_tran_err = maximum(skipmissing(tran_err))
			ok_tran = max_tran_err ≤ atol

		# 	BECAME_INFECTED logic
			ok_became = all(df.BECAME_INFECTED .== ((df.ACTIVE .== 1) .& (df.TRANS .== 1)))

		return (
			ok_prob_act = ok_prob,
			ok_tranprob = ok_tran,
			ok_became   = ok_became,
			max_abs_prob_err     = max_prob_err,
			max_abs_tranprob_err = max_tran_err
		)
	end

	β = 0.02
	checks = trace_sanity_checks(df; beta=β)
	println(checks)
	@assert checks.ok_prob_act "PROB_ACT != logistic(LWACT)"
	@assert checks.ok_tranprob "TRANPROB != 1 - (1-β)^EDGWGT"
	@assert checks.ok_became "BECAME_INFECTED != (ACTIVE==1 && TRANS==1)"

#	Attempt Coverage Check
	function trace_attempt_coverage(df::DataFrame)::Bool
		"""
		Args:
			df::DataFrame: Trace data
		Returns:
			Bool: true if within each (IT,T,J) group, maximum K equals NBRS exactly
		Notes:
			Assures we logged every neighbor attempt per ego/day.
		"""
		#	Group by (IT,T,J).
			g = groupby(df, [:IT,:T,:J])
		#	Check K coverage equals NBRS.
			all_ok = true
			for sub in g
				maxK = maximum(sub.K)
				nb   = first(sub.NBRS)
				if maxK != nb
					all_ok = false
					break
				end
			end
			return all_ok
	end

	coverage_ok = trace_attempt_coverage(df)
	println(("attempt_coverage" => coverage_ok))
	@assert coverage_ok "Trace is missing some neighbor attempts (K != NBRS)."

#	Build Fast Trace Index
	function build_trace_index(df::DataFrame)
		"""
		Args:
			df::DataFrame: Trace data
		Returns:
			Dict{NTuple{4,Int},NamedTuple}: key=(it,t,ego,colpos) → values for override
		Notes:
			Use COLPOS (not K) to align with SAS column order in s_alst[ego,].
		"""
		#	Select only fields needed at lookup sites.
			idx = Dict{NTuple{4,Int},NamedTuple}()
			for row in eachrow(df)
				key = (row.IT, row.T, row.EGO, row.COLPOS)
				val = (active = row.ACTIVE,
				       trans  = row.TRANS,
				       issymp = row.ISSYMPTOMATIC,
				       edgwgt = row.EDGWGT,
				       peers  = row.PEERSINF,
				       ninf   = row.NINF)
				idx[key] = val
			end
			return idx
	end

	idx = build_trace_index(df)
	length(idx)

#	Lookup Attempt From Index
	function trace_lookup(idx::Dict{NTuple{4,Int},NamedTuple}, it::Int, t::Int, ego::Int, colpos::Int)
		"""
		Args:
			idx::Dict: From `build_trace_index`
			it::Int: Iteration
			t::Int: Day
			ego::Int: Ego node id (SAS EGO)
			colpos::Int: Column position in SAS s_alst[ego,]
		Returns:
			NamedTuple: (active::Int, trans::Int, issymp::Int, edgwgt::Float64, peers::Int, ninf::Int)
		Notes:
			Throws a KeyError if the attempt isn’t in the trace.
		"""
		#	Strict retrieval.
			return idx[(it,t,ego,colpos)]
	end
                
#	Attempts For One Ego/Day
	function attempts_for_ego(df::DataFrame, it::Int, t::Int, ego::Int)::DataFrame
		"""
		Args:
			df::DataFrame: Trace data
			it,t,ego: Selection keys
		Returns:
			DataFrame: Attempts sorted by K (scan order)
		Notes:
			Useful for quick spot checks.
		"""
		#	Filter and sort by K.
			sub = filter(r -> r.IT==it && r.T==t && r.EGO==ego, df)
			sort!(sub, [:K])
			return sub
	end

	sub_584 = attempts_for_ego(df, 1, 1, 584)  # returns attempts sorted by K
	first(sub_584, 10)  

	if nrow(sub_584) > 0
		colpos1 = sub_584.COLPOS[1]
		trace_lookup(idx, 1, 1, 584, colpos1)
	end

#	Recompute LWACT/PROB_ACT For Audit
	function recompute_scores(alst::Matrix{Int64}, df::DataFrame;
		b_int::Float64,
		b_close::Float64,
		b_cxn_peers::Float64,
		b_cxn_total::Float64,
		b_cxn_symp::Float64,
		b_cls_x_smp::Float64)::NamedTuple
		"""
		Args:
			df::DataFrame: Trace data
			…: Model coefficients
		Returns:
			NamedTuple: (max_abs_lwact_diff, max_abs_prob_diff)
		Notes:
			Uses PEERSINF and NINF columns straight from the trace.
		"""
		#	Set Parameters
			n_const = length(unique(alst[:,1]))
			global_scale = Float64.(df.NINF) ./ (n_const / 3)

		#	Recompute.
			lw = b_int .+
				 b_cxn_symp .* Float64.(df.ISSYMPTOMATIC) .+
				 b_close     .* df.EDGWGT .+
				 b_cxn_peers .* Float64.(df.PEERSINF) .+
				 b_cxn_total .* global_scale .+
				 b_cls_x_smp .* Float64.(df.PEERSINF) .* Float64.(df.ISSYMPTOMATIC)

			prob = 1.0 ./ (1.0 .+ exp.(-lw))

			max_lw = maximum(abs.(lw .- df.LWACT))
			max_pb = maximum(abs.(prob .- df.PROB_ACT))

			return (max_abs_lwact_diff = max_lw,
			        max_abs_prob_diff  = max_pb)
	end

	b_int         = -0.1
	b_close       =  1.0
	b_cxn_peers   = -1.0
	b_cxn_total   = -6.50
	b_cxn_symp    = -1.5
	b_cls_x_smp   = -0.1
	audit = recompute_scores(alst, df;
		b_int=b_int, b_close=b_close, b_cxn_peers=b_cxn_peers,
		b_cxn_total=b_cxn_total, b_cxn_symp=b_cxn_symp, b_cls_x_smp=b_cls_x_smp)
	println(audit)

	by_day = combine(groupby(df, :T),
		nrow                 => :n_attempts,
		:ACTIVE             => sum => :n_activated,
		:BECAME_INFECTED    => sum => :n_transmitted,
	)

	by_it_day = combine(groupby(df, [:IT, :T]),
		nrow                 => :n_attempts,
		:ACTIVE             => sum => :n_activated,
		:BECAME_INFECTED    => sum => :n_transmitted,
	)

	transform!(by_it_day,
		[:n_activated, :n_attempts]   => ((a,n)->a./n) => :activation_rate,
		[:n_transmitted, :n_attempts] => ((x,n)->x./n) => :transmission_rate,
	)

# 	NINF should be constant within (IT,T)
	const_ninf = combine(groupby(df, [:IT, :T]), :NINF => (x -> length(unique(x))) => :n_unique_ninf)
	all(const_ninf.n_unique_ninf .== 1)

#	# of egos processed that day == NINF_t0
	egos_per_day = combine(groupby(df, [:IT, :T])) do sdf
    	(; egos_with_attempts = length(unique(sdf.EGO)), ninf_t0 = first(sdf.NINF))
	end

	diff = egos_per_day.ninf_t0 .- egos_per_day.egos_with_attempts
	println((
		all_nonnegative = all(diff .>= 0),      # should be true
		any_missing     = any(diff .> 0),       # true if some egos had no attempts
		summary         = describe(diff)
	))
	
	by_day = combine(groupby(df, [:IT,:T]),
		:EGO => (x->length(unique(x))) => :egos_with_attempts,
		:NINF => first => :ninf_t0,
	)
	by_day.diff = by_day.ninf_t0 .- by_day.egos_with_attempts

	println((
		all_nonnegative = all(by_day.diff .>= 0),
		any_silent_days = any(by_day.diff .> 0),
		n_silent_days   = count(>(0), by_day.diff),
	))

# 	Reconstruct I_t0 for each (IT,T) purely from the trace’s transmissions.
# 	Seed is taken as the first EGO seen on the first day for that iteration.
	function reconstruct_infected_sets(df::DataFrame; rec_t::Int=14)
		sets = Dict{Tuple{Int,Int}, Set{Int}}()

		for it in unique(df.IT)
			sdf = df[df.IT .== it, :]
			tmin = minimum(sdf.T)
			tmax = maximum(sdf.T)

			#	Seed guess: first EGO appearing on the earliest day
				seed = first(sdf.EGO[sdf.T .== tmin])

			#	Day of infection for everyone (seed is "pre-day", so treat as tmin-1)
				tinf = Dict{Int,Int}(seed => tmin - 1)
				for row in eachrow(sdf)
					if row.BECAME_INFECTED == 1
						# if someone gets infected multiple times in trace (shouldn’t), keep earliest
						tinf[row.ALTER] = get(tinf, row.ALTER, row.T)
					end
				end

			#	Build start-of-day infected sets
				for t in tmin:tmax
					S = Set{Int}()
					for (node, t0) in tinf
						# Node is infected at the start of day t if t ∈ (t0+1) … (t0+rec_t)
						if (t >= t0 + 1) && (t <= t0 + rec_t)
							push!(S, node)
						end
					end
					sets[(it,t)] = S
				end
		end

		return sets
	end

	sets = reconstruct_infected_sets(df; rec_t=14)

	cmp = DataFrame(IT = Int[], T = Int[], n_recon = Int[], ninf_t0 = Int[], diff = Int[])
	for it in unique(df.IT), t in unique(df.T[df.IT .== it])
		n_recon = length(sets[(it,t)])
		ninf_t0 = first(unique(df.NINF[(df.IT .== it) .&& (df.T .== t)]))
		push!(cmp, (it, t, n_recon, ninf_t0, n_recon - ninf_t0))
	end

	println((
		all_match = all(cmp.diff .== 0),
		n_mismatched_days = count(!=(0), cmp.diff),
	))

# 	Check monotone (non-decreasing) PEERSINF across K within (IT,T,EGO),
# 	and whether it steps up immediately after a success.
	function audit_within_day_peers(df::DataFrame)
		g = groupby(df, [:IT, :T, :EGO])
		n_groups = length(g)
		monotone_ok = 0
		step_ups = 0
		step_ups_follow_success = 0
		total_successes = 0

		for sdf in g
			#	Defining Parameters
				sort!(sdf, [:K])  # ensure inner-loop order
				p = sdf.PEERSINF
				b = sdf.BECAME_INFECTED

			# 	monotone non-decreasing across K?
				if all(Base.diff(p) .>= 0)
					monotone_ok += 1
				end

				d = Base.diff(p)
				step_ups += count(>(0), d)

				succ_idx = findall(==(1), b)
				total_successes += length(succ_idx)
				for si in succ_idx
					if si < nrow(sdf) && p[si+1] > p[si]
						step_ups_follow_success += 1
					end
				end
		end

		return (
			groups = n_groups,
			groups_monotone = monotone_ok,
			share_monotone = n_groups == 0 ? NaN : monotone_ok / n_groups,
			step_ups = step_ups,
			total_successes = total_successes,
			step_ups_follow_success = step_ups_follow_success,
		)
	end

	audit_within_day_peers(df)

#	Check 1: PEERSINF persistence across days
	function check_peersinf_persistence(df::DataFrame)
		# 	Group by IT and ALTER
			grouped = groupby(df, [:IT, :ALTER])
		
			alters_multi_day = 0
			any_decrease = false
			max_peers = 0
		
		#	Looping over the Grouped Dataset
			for grp in grouped
				#	Get unique days for this group
					days = sort(unique(grp.T))
				
					if length(days) > 1
						alters_multi_day += 1
						
						# Get first PEERSINF value for each day
						first_vals = Int[]
						for d in days
							day_rows = grp[grp.T .== d, :]
							if nrow(day_rows) > 0
								push!(first_vals, day_rows.PEERSINF[1])
							end
						end
						
						# Check for decreases
						if length(first_vals) > 1
							for i in 2:length(first_vals)
								if first_vals[i] < first_vals[i-1]
									any_decrease = true
									break
								end
							end
						end
					end
				
				# 	Track max PEERSINF
					if length(grp.PEERSINF) > 0
						max_peers = max(max_peers, maximum(grp.PEERSINF))
					end
			end
		
		#	Return Results
			return (
				alters_seen_multiple_days = alters_multi_day,
				any_alter_sees_decrease = any_decrease,
				max_peersinf_observed = max_peers
			)
	end
	persistence_check = check_peersinf_persistence(df)
	println("PEERSINF Persistence Check:")
	println(persistence_check)

#	Check 2: Verify ego's neighbors get updated on transmission
	function check_neighbor_updates(df::DataFrame)
		"""
		When ego successfully transmits, do ALL ego's future attempts 
		show increased PEERSINF?
		"""
		transmissions = df[df.BECAME_INFECTED .== 1, :]
		
		updates_found = 0
		updates_expected = 0
		
		for row in eachrow(transmissions[1:min(10, nrow(transmissions)), :])
			# Find ego's subsequent attempts same day
			post_trans = df[(df.IT .== row.IT) .& 
					       (df.T .== row.T) .& 
					       (df.EGO .== row.EGO) .& 
					       (df.K .> row.K), :]
			
			if nrow(post_trans) > 0
				# Did PEERSINF increase for remaining neighbors?
				increases = post_trans.PEERSINF .> row.PEERSINF
				updates_found += sum(increases)
				updates_expected += nrow(post_trans)
			end
		end
		
		return (
			updates_found = updates_found,
			updates_expected = updates_expected,
			update_rate = updates_expected > 0 ? updates_found/updates_expected : NaN
		)
	end

	neighbor_update_check = check_neighbor_updates(df)
	println("\nNeighbor Update Check:")
	println(neighbor_update_check)

# 	Find a transmission event and trace what happens
	function trace_transmission_example(df::DataFrame)
		# 	Find first transmission
			trans = df[df.BECAME_INFECTED .== 1, :]
			if nrow(trans) == 0
				return nothing
			end
		
		# 	Take first transmission
			t_row = trans[1, :]
		
			println("Transmission Example:")
			println("  IT=$(t_row.IT), T=$(t_row.T), EGO=$(t_row.EGO) -> ALTER=$(t_row.ALTER)")
			println("  At K=$(t_row.K) with PEERSINF=$(t_row.PEERSINF)")
			
		# 	Get all of ego's attempts that day
			ego_day = df[(df.IT .== t_row.IT) .& (df.T .== t_row.T) .& (df.EGO .== t_row.EGO), :]
			sort!(ego_day, :K)
			
			println("\nEgo's attempts that day:")
			for row in eachrow(ego_day)
				marker = row.K == t_row.K ? " <- TRANSMISSION HERE" : ""
				println("  K=$(row.K): ALTER=$(row.ALTER), PEERSINF=$(row.PEERSINF)$marker")
			end
		
		#	Return Ego Trace
			return ego_day
	end

	example = trace_transmission_example(df)

# 	Check 3: Verify that infected nodes stop appearing as EGO
	function check_infected_stop_transmitting(df::DataFrame)
		"""
		Once someone recovers (after 14 days), do they stop appearing as EGO?
		"""
		# 	Track when each node was infected and when they last appeared as ego
			ego_history = combine(groupby(df, [:IT, :EGO])) do sdf
				first_day = minimum(sdf.T)
				last_day = maximum(sdf.T)
				got_infected = any(sdf[sdf.T .== first_day, :].ALTER .∈ 
								Ref(df[(df.IT .== first(sdf.IT)) .& 
										(df.BECAME_INFECTED .== 1), :ALTER]))
				DataFrame(
					first_ego_day = first_day,
					last_ego_day = last_day,
					ego_duration = last_day - first_day + 1,
					was_seed = first_day == 1
				)
			end
		
		# 	Check if ego duration matches recovery time
			expected_duration = 14  # recovery time
			durations = ego_history[.!ego_history.was_seed, :ego_duration]
			
		return (
			mean_ego_duration = mean(durations),
			max_ego_duration = maximum(durations),
			all_stop_by_14 = all(durations .<= expected_duration)
		)
	end

# 	Check 4: Verify global effect calculation
	function check_global_effect(df::DataFrame, n_nodes::Int)
		"""
		Verify that NINF / (n/3) matches our understanding
		"""
		# 	Check unique NINF values per day
			ninf_by_day = combine(groupby(df, [:IT, :T]), 
								:NINF => unique => :ninf_values)
		
		# 	Calculate what global effect should be
			expected_scale = n_nodes / 3
		
		# 	Sample some calculations
			sample = df[1:min(100, nrow(df)), :]
			sample.global_effect = sample.NINF ./ expected_scale
		
			return (
				all_ninf_consistent = all(length.(ninf_by_day.ninf_values) .== 1),
				expected_divisor = expected_scale,
				sample_global_effects = first(sample.global_effect, 10)
			)
	end

# 	Check 5: Verify infection order processing
	function check_infection_order(df::DataFrame)
		"""
		Are infected nodes processed in a consistent order each day?
		"""
		# For each (IT, T), get the order of EGOs
		ego_order = combine(groupby(df, [:IT, :T])) do sdf
			egos = Int[]
			current_ego = 0
			for row in eachrow(sort(sdf, [:J, :K]))
				if row.EGO != current_ego
					push!(egos, row.EGO)
					current_ego = row.EGO
				end
			end
			DataFrame(ego_sequence = [egos])
		end
		
		# Check if order is consistent (e.g., always ascending by ID)
		ascending_ordered = 0
		for row in eachrow(ego_order)
			if length(row.ego_sequence) > 1
				if issorted(row.ego_sequence)
					ascending_ordered += 1
				end
			end
		end
		
		return (
			total_days = nrow(ego_order),
			days_with_ascending_ego_order = ascending_ordered,
			sample_ego_sequence = ego_order[1:min(5, nrow(ego_order)), :]
		)
	end

# 	Run Checks 3-5
	println("Infected Stop Transmitting Check:")
	println(check_infected_stop_transmitting(df))

	println("\nGlobal Effect Check (n_nodes = 614):")
	println(check_global_effect(df, 614))

	println("\nInfection Order Check:")
	println(check_infection_order(df))

####################
#   PACKAGE TESTS  #
####################

#	Setting Working Directory to the Test Data Directory
	cd("/workspace/data/sim_test_data/")
	pwd()

#	Test 1: Help (single & sweep)
	diffusion_sim.CLI.cli_main(["single", "--help"])
	diffusion_sim.CLI.cli_main(["sweep",  "--help"])
	diffusion_sim.CLI.cli_main(["multi",  "--help"])

#	Test 2: Single run end-to-end
	single_out = "/workspace/data/sim_test_data/result_single.json"

	single_args = [
		"single",
		"--graphml", "/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml",
		# "--sas-transformation",            # uncomment if you want SAS preprocessing
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
		"--out", single_out
	]
	diffusion_sim.CLI.cli_main(single_args)

	single_csv = replace(single_out, r"\.json$" => "") * "_infection_log.csv"
	@assert isfile(single_csv) "Single-run CSV not found at $(single_csv)"
	df_single = CSV.read(single_csv, DataFrame)
	println("single rows: ", nrow(df_single), "  cols: ", names(df_single))
	df_single[1:10,:]

# 	Test 3: Sweep run (comma-separated vectors)
	sweep_out = "/workspace/ergm-based-diffusion-simulations/result_sweep.json"

	sweep_args = [
		"sweep",
		"--graphml", "/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml",
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
		"--out", sweep_out
	]
	diffusion_sim.CLI.cli_main(sweep_args)

	sweep_csv = replace(sweep_out, r"\.json$" => "") * "_infection_log.csv"
	@assert isfile(sweep_csv) "Sweep CSV not found at $(sweep_csv)"
	df_sweep = CSV.read(sweep_csv, DataFrame)
	println("sweep rows: ", nrow(df_sweep), "  cols: ", names(df_sweep))

# 	Test 4: Sweep columns & grouping sanity
	expected = [:time, :n_infected, :prop_ever_infected, :prop_currently_infected, :n_recovered,
 				:prop_recovered, :b_cxn_peers, :b_cxn_total,:b_cxn_symp,:iter]
	present = names(df_single)
	print(expected)
	print(present)

# 	Each (bp,bg,bs,iter) should produce a full trajectory
	by_keys = groupby(df_sweep, [:b_cxn_peers, :b_cxn_total, :b_cxn_symp, :iter])
	traj_lengths = combine(by_keys, :time => length => :nrows)
	println("distinct trajectories: ", nrow(traj_lengths))
	println(first(traj_lengths, 8))

# 	Check each group’s min time == 0
	g = groupby(df_sweep, [:b_cxn_peers, :b_cxn_total, :b_cxn_symp, :iter])
	mins = combine(g, :time => minimum => :tmin)
	@assert all(mins.tmin .== 0.0) "Some trajectories do not start at time==0"

# 	Check that no time exceeds your horizon (e.g., 200)
	@assert maximum(df_sweep.time) <= 200 "Found time > maxtime"

# 	Test 5: Processsing 10 GraphML files from the directory
	dir = "/workspace/data/sim_test_data/synthetic_village_m_graphml"
	all_files = readdir(dir, join=true)
	graphml_files = filter(f -> endswith(f, ".graphml"), all_files)[1:10]

	file_list = join(graphml_files, ",")

	multi_out = "/workspace/ergm-based-diffusion-simulations/multi_test_results.json"
	args = [
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
		"--out", multi_out
	]

	diffusion_sim.CLI.cli_main(args)

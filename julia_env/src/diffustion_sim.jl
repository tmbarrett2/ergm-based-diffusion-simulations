module diffustion_sim
#   Packages
    using CSV
    using DataFrames
	using EzXML
	using Random
    using RCall
    using Statistics
	using StatsBase
   
#   SUPPORT FUNCTIONS

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

#	Helper Function for sim_prep: load and parse GraphML file
	function load_graphml(filepath::String)
		"""
		Args:
			filepath::String: path to .graphml file
		Returns:
			NamedTuple: (nodes, edges, weights, node_labels) extracted from GraphML
		Notes:
			Parses GraphML directly to extract network structure.
			Handles both weighted and unweighted networks.
		"""
		
		#	Read and parse XML
			doc = readxml(filepath)
			root = doc.root
			
		#	Find graph element
			graph = nothing
			for child in eachelement(root)
				if child.name == "graph"
					graph = child
					break
				end
			end
			
			if graph === nothing
				error("No graph element found in GraphML file")
			end
			
		#	Determine if directed
			is_directed = haskey(graph, "edgedefault") && graph["edgedefault"] == "undirected" ? false : true
			
		#	Extract nodes by direct traversal
			nodes_xml = []
			node_ids = String[]
			node_mapping = Dict{String, Int}()
			
			for elem in eachelement(graph)
				if elem.name == "node"
					push!(nodes_xml, elem)
					node_id = elem["id"]
					push!(node_ids, node_id)
					node_mapping[node_id] = length(node_ids)
				end
			end
			
		#	Extract edges by direct traversal
			edges_xml = []
			el1 = Int[]
			el2 = Int[]
			elwgt = Float64[]
			
			for elem in eachelement(graph)
				if elem.name == "edge"
					source = elem["source"]
					target = elem["target"]
					
					#	Map string IDs to integer indices
						if haskey(node_mapping, source) && haskey(node_mapping, target)
							src_idx = node_mapping[source]
							tgt_idx = node_mapping[target]
							
							push!(el1, src_idx)
							push!(el2, tgt_idx)
							
							#	Look for weight in data elements
								weight = 1.0
								for data_elem in eachelement(elem)
									if data_elem.name == "data" && haskey(data_elem, "key")
										# Check for various weight key names
										if data_elem["key"] in ["e_weight", "weight", "value", "d1"]
											weight_str = nodecontent(data_elem)
											if !isempty(weight_str)
												weight = parse(Float64, weight_str)
											end
											break
										end
									end
								end
								push!(elwgt, weight)
							
							#	Add reverse edge if undirected
								if !is_directed
									push!(el1, tgt_idx)
									push!(el2, src_idx)
									push!(elwgt, weight)
								end
						else
							println("Warning: Edge references non-existent node(s): $source -> $target")
						end
				end
			end
			
		#	Clean up
			EzXML.finalize(doc) 
			
		#	Check if we found data
			if isempty(node_ids)
				error("No nodes found in GraphML file")
			end
			if isempty(el1)
				error("No edges found in GraphML file")
			end
			println("Loaded network: $(length(node_ids)) nodes, $(length(el1)) directed edges")
			
		#	Return extracted data
			return (
				nodes = collect(1:length(node_ids)),
				edges = (el1, el2),
				weights = elwgt,
				node_labels = node_ids
			)
	end

#	Helper Function for sim_prep: convert valued edgelist to adjacency matrix
	function el2adjval(valid::Vector{Int}, el1::Vector{Int}, el2::Vector{Int}, elwgt::Vector{Float64})
		"""
		Args:
			valid::Vector{Int}: unique node identifiers
			el1::Vector{Int}: source nodes of edges
			el2::Vector{Int}: target nodes of edges
			elwgt::Vector{Float64}: edge weights
		Returns:
			Tuple{Matrix{Float64}, Vector{Int}}: weighted adjacency matrix and node ordering
		Notes:
			Converts edge list representation to matrix form.
			Matrix rows/columns ordered by valid node set.
		"""
		
		#	Initialize adjacency matrix
			nodeset = unique(valid)
			n = length(nodeset)
			adjmat = zeros(Float64, n, n)
			
		#	Create node index mapping
			node_to_idx = Dict(node => idx for (idx, node) in enumerate(nodeset))
		
		#	Populate adjacency matrix
			for i in 1:length(el1)
				iv = el1[i]
				jv = el2[i]
				wv = elwgt[i]
				
				#	Map nodes to matrix indices
					if haskey(node_to_idx, iv) && haskey(node_to_idx, jv)
						iloc = node_to_idx[iv]
						jloc = node_to_idx[jv]
						adjmat[iloc, jloc] = wv
					end
			end
		
		#	Return matrix with node ordering
			return (adjmat, nodeset)
	end

#	Helper Function for preprocess_network_sas_style: convert adjacency to matrix
	function adjacency_list_to_matrix(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: adjacency list (node ID in column 1, neighbors in 2:end)
			vlst::Matrix{Float64}: value list (matching structure)
		Returns:
			Tuple{Matrix{Float64}, Vector{Int}}: (adjacency_matrix, node_ids)
		Notes:
			Converts adjacency list format to square adjacency matrix.
		"""
		
		#	Get node IDs from first column
			node_ids = alst[:, 1]
			n = length(node_ids)
			
		#	Initialize adjacency matrix
			adj_matrix = zeros(Float64, n, n)
			
		#	Fill matrix from adjacency list
			for i in 1:n
				for j in 2:size(alst, 2)
					neighbor = alst[i, j]
					if neighbor > 0
						#	Find neighbor's row index
							neighbor_idx = findfirst(x -> x == neighbor, node_ids)
							if neighbor_idx !== nothing
								adj_matrix[i, neighbor_idx] = vlst[i, j]
							end
					end
				end
			end
			
		#	Return matrix and node ordering
			return (adj_matrix, node_ids)
	end

#	Helper Function for preprocess_network_sas_style: convert matrix back to adjacency list
	function matrix_to_adjacency_list(adj_matrix::Matrix{Float64}, node_ids::Vector{Int})
		"""
		Args:
			adj_matrix::Matrix{Float64}: weighted adjacency matrix
			node_ids::Vector{Int}: node IDs corresponding to matrix rows/columns
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (alst, vlst)
		Notes:
			Converts matrix format back to adjacency list format.
		"""
		
		#	Calculate maximum degree
			n = size(adj_matrix, 1)
			max_degree = maximum(sum(adj_matrix .> 0, dims=2))
			
		#	Initialize adjacency and value lists
			alst = zeros(Int, n, Int(max_degree) + 1)
			vlst = zeros(Float64, n, Int(max_degree) + 1)
			
		#	Fill node IDs in first column
			alst[:, 1] = node_ids
			vlst[:, 1] = Float64.(node_ids)
			
		#	Fill adjacency and value lists
			for i in 1:n
				neighbors = findall(x -> x > 0, adj_matrix[i, :])
				if !isempty(neighbors)
					for (col, j) in enumerate(neighbors)
						if col + 1 <= size(alst, 2)
							alst[i, col + 1] = node_ids[j]
							vlst[i, col + 1] = adj_matrix[i, j]
						end
					end
				end
			end
			
		#	Return lists
			return (alst, vlst)
	end

#	Preprocess Network to Match SAS IML Code
	function preprocess_network_sas_style(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: raw adjacency list from sim_prep
			vlst::Matrix{Float64}: raw value list from sim_prep  
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (processed_alst, processed_vlst)
		Notes:
			Replicates SAS preprocessing:
			1. Symmetrizes the network
			2. Adds common third-party connections
			3. Takes square root of common thirds
			4. Returns in adjacency list format
		"""
		
		#	Convert to matrix format
			adj_matrix, node_ids = adjacency_list_to_matrix(alst, vlst)
			
		#	Symmetrize: symmat = (amat + amat`)
			symmat = adj_matrix + transpose(adj_matrix)
			
		#	Calculate common thirds: com3rds = ((symmat>0)*(symmat>0))#(symmat>0)
		#	This is the element-wise product of the binary adjacency matrix with itself
		#	In SAS: A#B is element-wise multiplication, * is matrix multiplication
			binary_mat = Float64.(symmat .> 0)
			com3rds = (binary_mat * binary_mat) .* binary_mat
			
		#	Take square root of common thirds: com3rds = com3rds##(0.5)
		#	In SAS: ## is element-wise power
			com3rds = sqrt.(com3rds)
			
		#	Add common thirds to symmetrized matrix
			symmat = symmat + com3rds
			
		#	Convert back to adjacency list format
			processed_alst, processed_vlst = matrix_to_adjacency_list(symmat, node_ids)
			
		#	Return processed networks
			return (processed_alst, processed_vlst)
	end
	@doc raw"""
	**Description**
	Preprocesses a network to match the exact transformations performed in the SAS IML code.
	This ensures Julia simulations use the same network structure as SAS simulations.

	**Usage**
	`preprocess_network_sas_style(alst, vlst)`

	**Arguments**
	- `alst::Matrix{Int}`: Raw adjacency list from `sim_prep`
	- `vlst::Matrix{Float64}`: Raw value list from `sim_prep`

	**Details**
	The function performs the following transformations to match SAS:
	
	1. **Symmetrization**: Creates an undirected network by adding the adjacency matrix 
	   to its transpose: `symmat = A + A'`
	
	2. **Common third-party connections**: Identifies paths of length 2 between nodes.
	   For nodes i and j, this counts how many nodes k exist where both i→k and k→j.
	   Computed as: `com3rds = (A_binary * A_binary) ⊙ A_binary`
	
	3. **Square root weighting**: Takes element-wise square root of common thirds
	   to moderate their influence: `com3rds = sqrt(com3rds)`
	
	4. **Final network**: Adds common thirds to symmetrized network:
	   `final = symmat + com3rds`
	
	This preprocessing strengthens connections between nodes that share common neighbors,
	which can significantly affect epidemic dynamics by creating additional transmission
	pathways and stronger clustering.

	**Value**
	Returns a tuple `(processed_alst, processed_vlst)`:
	- `processed_alst::Matrix{Int}`: Preprocessed adjacency list
	- `processed_vlst::Matrix{Float64}`: Preprocessed edge weights

	**Examples**
	```julia
	# Load raw network
	network_data = sim_prep("network.graphml")
	
	# Apply SAS-style preprocessing
	alst_processed, vlst_processed = preprocess_network_sas_style(
	    network_data.alst, 
	    network_data.vlst
	)
	
	# Use in simulation (now matches SAS)
	results = sirdif(alst_processed, vlst_processed, [1], 
	                0.02, 14, 200, 0.75,
	                -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
	
	# Compare network properties
	raw_edges = sum(network_data.vlst .> 0)
	processed_edges = sum(vlst_processed .> 0)
	println("Raw edges: $raw_edges, Processed edges: $processed_edges")
	```

	**See Also**
	`sim_prep`, `sirdif`, `replicate_sas_simulation`
	""" preprocess_network_sas_style

#	SIMULATION FUNCTIONS

#	Network Data Preparation for SIR Simulation from GraphML
	function sim_prep(graphml_file::String; sas_transformation::Bool = false)
		"""
		Args:
			graphml_file::String: path to .graphml file
			sas_transformation::Bool: apply SAS-style preprocessing (default false)
		Returns:
			NamedTuple: (alst=adjacency_list, vlst=value_list)
		Notes:
			Loads GraphML file and converts to adjacency list format for sir_diffusion.
			If sas_transformation=true, applies symmetrization and common thirds.
		"""
		
		#	Load network from GraphML
			network_data = load_graphml(graphml_file)
			nl = network_data.nodes
			el1, el2 = network_data.edges
			elwgt = network_data.weights
			n = length(nl)
		
		#	Apply SAS transformation if requested
			if sas_transformation
				#	CRITICAL: load_graphml already doubles edges for undirected graphs
				#	We need to deduplicate to match what SAS pajread produces
				#	Keep only unique undirected edges
					edge_set = Set{Tuple{Int,Int}}()
					el1_dedup = Int[]
					el2_dedup = Int[]
					
					for (src, tgt) in zip(el1, el2)
						edge_key = src <= tgt ? (src, tgt) : (tgt, src)
						if !(edge_key in edge_set)
							push!(edge_set, edge_key)
							push!(el1_dedup, src)
							push!(el2_dedup, tgt)
						end
					end
				
				#	Build adjacency matrix exactly as SAS adj() does
					unique_nodes = sort(unique(vcat(el1_dedup, el2_dedup)))
					n = length(unique_nodes)
					node_to_idx = Dict(node => idx for (idx, node) in enumerate(unique_nodes))
					
				#	Create adjacency matrix (counting edges)
					amat = zeros(Float64, n, n)
					for (sender, receiver) in zip(el1_dedup, el2_dedup)
						i = node_to_idx[sender]
						j = node_to_idx[receiver]
						amat[i, j] += 1.0  # SAS adj() increments for each edge
					end
				
				#	Step 2a: Symmetrize (amat + amat')
					symmat = amat + transpose(amat)
					
				#	Step 2b: Calculate common thirds
				#	com3rds = ((symmat>0)*(symmat>0))#(symmat>0)
					binary = Float64.(symmat .> 0)
					com3rds = (binary * binary) .* binary
					
				#	Step 2c: Square root of common thirds
				#	com3rds = com3rds##0.5
					com3rds = sqrt.(com3rds)
					
				#	Step 2d: Add to symmetrized matrix
				#	symmat = symmat + com3rds
					symmat = symmat + com3rds
					
				#	Convert to adjacency list format
					max_degree = Int(maximum(sum(symmat .> 0, dims=2)))
					alst = zeros(Int, n, max_degree + 1)
					vlst = zeros(Float64, n, max_degree + 1)
					
				#	Fill ID columns (using the sorted unique nodes)
					alst[:, 1] = unique_nodes
					vlst[:, 1] = Float64.(unique_nodes)
					
				#	Fill neighbor lists
					for i in 1:n
						neighbors = findall(x -> x > 0, symmat[i, :])
						for (col_idx, j) in enumerate(neighbors)
							if col_idx + 1 <= size(alst, 2)
								alst[i, col_idx + 1] = unique_nodes[j]
								vlst[i, col_idx + 1] = symmat[i, j]
							end
						end
					end
					
					println("Applied SAS transformation: symmetrization + common thirds")
			else
				#	No SAS transformation - use raw network
				#	Create binary adjacency matrix for structure
					adj_mat_binary, nodeset = el2adjval(nl, el1, el2, ones(Float64, length(elwgt)))
					
				#	Create weighted adjacency matrix
					adj_mat_weighted, _ = el2adjval(nl, el1, el2, elwgt)
					
				#	Calculate maximum degree
					max_degree = Int(maximum(sum(adj_mat_binary .> 0, dims=2)))
					
				#	Initialize adjacency list and value list
					alst = zeros(Int, n, max_degree + 1)  # +1 for node ID column
					vlst = zeros(Float64, n, max_degree + 1)
					
				#	Populate node IDs in first column
					for i in 1:n
						alst[i, 1] = nodeset[i]
						vlst[i, 1] = Float64(nodeset[i])
					end
				
				#	Fill adjacency and value lists
					for i in 1:n
						#	Find neighbors
							neighbors = findall(x -> x > 0, adj_mat_binary[i, :])
							
							if !isempty(neighbors)
								#	Map indices to node IDs and extract weights
									neighbor_ids = [nodeset[j] for j in neighbors]
									weights = [adj_mat_weighted[i, j] for j in neighbors]
									
								#	Fill adjacency list
									for (col, nid) in enumerate(neighbor_ids)
										if col + 1 <= size(alst, 2)
											alst[i, col + 1] = nid
										end
									end
									
								#	Fill value list
									for (col, wgt) in enumerate(weights)
										if col + 1 <= size(vlst, 2)
											vlst[i, col + 1] = wgt
										end
									end
							end
					end
			end
		
		#	Return formatted lists for simulation
			return (alst = alst, vlst = vlst)
	end
	@doc raw"""
	**Description**
	Prepares network data from a GraphML file for use in SIR diffusion simulations.
	Directly parses GraphML format and converts to adjacency list representation with
	optional SAS-style preprocessing for comparison with SAS PROC IML simulations.

	**Usage**
	`sim_prep(graphml_file; sas_transformation=false)`

	**Arguments**
	- `graphml_file::String`: Path to .graphml file containing network data
	- `sas_transformation::Bool`: Apply SAS-style preprocessing (default `false`)

	**Details**
	The function performs the following operations:
	1. Parses GraphML file using EzXML to extract nodes, edges, and weights
	2. Handles both directed and undirected graphs automatically
	3. If `sas_transformation=true`, applies these transformations:
	- **Symmetrization**: Creates undirected network via `A + A'`
	- **Common thirds**: Adds edges between nodes sharing neighbors via `(A*A) ⊙ A`
	- **Square root weighting**: Applies `sqrt()` to common third connections
	- **Combination**: Final network = symmetrized + sqrt(common thirds)
	4. Converts to adjacency list format where:
	- Each row represents a node
	- First column contains the node ID
	- Subsequent columns contain neighbor IDs (0 for empty slots)
	5. Creates parallel value list containing edge weights

	The SAS transformation significantly changes network topology by adding edges
	between nodes that share common neighbors, which can increase epidemic spread
	by creating additional transmission pathways.

	GraphML weight attributes are detected automatically. The function looks for
	edge attributes named "weight", "value", "e_weight", or "d1". If no weights 
	are found, all edges are assigned weight 1.0.

	**Value**
	Returns a NamedTuple with two fields:
	- `alst`: Matrix{Int} where row i contains node i's ID followed by its neighbor IDs
	- `vlst`: Matrix{Float64} with same structure containing edge weights

	**Examples**
	```julia
	# Load network without SAS preprocessing (raw network)
	network_data = sim_prep("network.graphml")

	# Load network with SAS preprocessing (for comparison with SAS)
	network_sas = sim_prep("network.graphml", sas_transformation=true)

	# Compare network properties
	raw_edges = sum(network_data.vlst .> 0)
	sas_edges = sum(network_sas.vlst .> 0)
	println("Raw edges: $raw_edges, SAS-processed edges: $sas_edges")

	# Use in simulations
	# For Julia-only studies:
	results_julia = sirdif(network_data.alst, network_data.vlst, [1], 
						0.02, 14, 200, 0.75, -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)

	# For SAS comparison studies:
	results_sas_style = sirdif(network_sas.alst, network_sas.vlst, [1], 
							0.02, 14, 200, 0.75, -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
	```

	**See Also**
	`sirdif`, `replicate_sas_simulation`, `sas_simulation_comparer`
	""" sim_prep

#	Helper Function for sirdif: convert matrix adjacency to vector format
	function matrix_to_vector_adjacency(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: adjacency matrix from sim_prep
			vlst::Matrix{Float64}: value matrix from sim_prep
		Returns:
			Tuple{Vector{Vector{Int}}, Vector{Vector{Float64}}}: (alst_vec, vlst_vec)
		Notes:
			Converts matrix format to vector format for efficiency.
			Pre-allocates output vectors.
		"""
		
		#	Get dimensions
			n = size(alst, 1)
			
		#	Pre-allocate vector versions
			alst_vec = Vector{Vector{Int}}(undef, n)
			vlst_vec = Vector{Vector{Float64}}(undef, n)
			
		#	Convert each row
			for i in 1:n
				#	Count non-zero neighbors first
					n_neighbors = count(x -> x > 0, @view alst[i, 2:end])
					
				#	Pre-allocate neighbor and weight vectors
					neighbors = Vector{Int}(undef, n_neighbors)
					weights = Vector{Float64}(undef, n_neighbors)
					
				#	Fill pre-allocated vectors
					idx = 1
					for j in 2:size(alst, 2)
						if alst[i, j] > 0
							neighbors[idx] = alst[i, j]
							weights[idx] = vlst[i, j]
							idx += 1
						end
					end
					
				#	Assign to output
					alst_vec[i] = neighbors
					vlst_vec[i] = weights
			end
			
		#	Return converted format
			return (alst_vec, vlst_vec)
	end

#	SIR Diffusion Simulation with Social Feedback
	function sirdif(
		alst::Union{Matrix{Int}, Vector{Vector{Int}}},
		vlst::Union{Matrix{Float64}, Vector{Vector{Float64}}},
		infectedp::Vector{Int}, inf_r::Float64, rec_t::Int,
		maxtime::Int, p_symp::Float64, b_int::Float64,
		b_close::Float64, b_cxn_peers::Float64, b_cxn_total::Float64,
		b_cxn_symp::Float64, b_cls_x_smp::Float64;
		transmission_method::Symbol = :weighted,
		immunity_duration::Union{Int,Nothing} = nothing
	)
		"""
		Args:
			alst: Adjacency list (Matrix or Vector{Vector})
			vlst: Edge weights aligned to alst
			infectedp: Initial infected node IDs
			inf_r: Base transmission probability
			rec_t: Recovery time in days
			maxtime: Maximum simulation days
			p_symp: Probability infected is symptomatic
			b_int: Baseline interaction coefficient
			b_close: Edge weight coefficient
			b_cxn_peers: Peer infection coefficient
			b_cxn_total: Global infection coefficient
			b_cxn_symp: Symptomatic coefficient
			b_cls_x_smp: Peers × symptomatic interaction
			transmission_method: :weighted or :simple (default :weighted)
			immunity_duration: Days of immunity after recovery (default nothing)
		Returns:
			Dict with infection_log, total_time, final_state
		Notes:
			Exact SAS behavior: FIFO ordering, cumulative peer counts, and single-row padding at extinction.
		"""

		#	Cache network size
			if alst isa Matrix{Int}
				n_nodes = size(alst, 1)
			else
				n_nodes = length(alst)
			end

		#	Normalize inputs to vector-of-vectors; collect node ids
			if alst isa Matrix{Int}
				alst_vec   = Vector{Vector{Int}}(undef, n_nodes)
				vlst_vec   = Vector{Vector{Float64}}(undef, n_nodes)
				unique_ids = Vector{Int}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					ego_id        = alst[i, 1]
					unique_ids[i] = ego_id
					row_ids       = alst[i, 2:end]
					row_wgts      = vlst[i, 2:end]
					nz_mask       = row_ids .> 0
					neighbors     = row_ids[nz_mask]
					weights       = row_wgts[nz_mask]
					alst_vec[i]   = [ego_id; neighbors]
					vlst_vec[i]   = [Float64(ego_id); weights]
				end
			else
				alst_vec   = Vector{Vector{Int}}(undef, n_nodes)
				vlst_vec   = Vector{Vector{Float64}}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					alst_vec[i] = alst[i]
					vlst_vec[i] = vlst[i]
				end
				unique_ids = Vector{Int}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					unique_ids[i] = alst_vec[i][1]
				end
			end

		#	Map node id → row index
			id_to_idx = Dict{Int,Int}(unique_ids[i] => i for i in 1:n_nodes)

		#	Start wall clock
			start_time = time()

		#	State matrix: [id, I, S, R, t_rec, nbrsinf, infection_order]
			state = zeros(Int, n_nodes, 7)
			@inbounds begin
				state[:, 1] .= unique_ids
				state[:, 3] .= 1
			end

		#	Column indices
			ID_COL               = 1
			INFECTED_COL         = 2
			SUSCEPTIBLE_COL      = 3
			RECOVERED_COL        = 4
			TIME_TO_RECOVERY_COL = 5
			NBRSINF_COL          = 6
			INFECTION_ORDER_COL  = 7

		#	Initialize infection order counter
			infection_counter = 0

		#	Seed infections
			@inbounds for inf_id in infectedp
				idx = id_to_idx[inf_id]
				state[idx, INFECTED_COL]         = 1
				state[idx, SUSCEPTIBLE_COL]      = 0
				state[idx, TIME_TO_RECOVERY_COL] = rec_t
				infection_counter += 1
				state[idx, INFECTION_ORDER_COL]  = infection_counter
			end

		#	Susceptible adjacency (preserve original column order)
			s_alst_vec = Vector{Vector{Int}}(undef, n_nodes)
			s_vlst_map = Vector{Dict{Int,Float64}}(undef, n_nodes)
			@inbounds for i in 1:n_nodes
				neighbors   = alst_vec[i][2:end]
				weights     = vlst_vec[i][2:end]
				s_alst_vec[i] = copy(neighbors)
				s_vlst_map[i] = Dict(zip(neighbors, weights))
			end

		#	In-place vector filter helper
			@inline function remove_id!(v::Vector{Int}, id::Int)
				w = 1
				@inbounds for i in 1:length(v)
					x = v[i]
					if x != id
						v[w] = x
						w += 1
					end
				end
				if w <= length(v)
					resize!(v, w - 1)
				end
				return nothing
			end

		#	Remove initially infected from all susceptible lists
			@inbounds for inf_id in infectedp
				for i in 1:n_nodes
					remove_id!(s_alst_vec[i], inf_id)
				end
			end

		#	Persistent peer-infected counts (cumulative)
			nbrsinf = zeros(Int, n_nodes)
			@inbounds for i in 1:n_nodes
				if state[i, INFECTED_COL] == 1
					for alter_id in s_alst_vec[i]
						alter_idx = id_to_idx[alter_id]
						nbrsinf[alter_idx] += 1
					end
				end
			end

		#	Time series buffer (preallocated; may finish early)
			timesum = Matrix{Float64}(undef, maxtime + 1, 7)
			wrow    = 1
			n_initial = length(infectedp)
			pinf      = n_initial / n_nodes
			timesum[wrow, :] = Float64[0.0, n_initial, pinf, 0.0, 0.0, 0.0, 0.0]

		#	Per-day buffers
			infected_idx_buf   = Vector{Int}(undef, n_nodes)
			infected_order_buf = Vector{Int}(undef, n_nodes)
			infected_perm_buf  = Vector{Int}(undef, 0)
			ordered_idx_buf    = Vector{Int}(undef, n_nodes)

		#	Main loop over days
			@inbounds for t in 1:maxtime

				#	Collect infected at start-of-day
					ninf = 0
					for i in 1:n_nodes
						if state[i, INFECTED_COL] == 1
							ninf += 1
							infected_idx_buf[ninf]   = i
							infected_order_buf[ninf] = state[i, INFECTION_ORDER_COL]
						end
					end

				#	Extinction padding (SAS: single row at maxtime)
					if ninf == 0
						NR_now = 0
						for i in 1:n_nodes
							NR_now += state[i, RECOVERED_COL]
						end
						NI_now = 0
						prop_ever_now = NR_now == 0 ? 0.0 : NR_now / n_nodes
						prop_cur_now  = 0.0
						prop_rec_now  = NR_now == 0 ? 0.0 : 1.0
						wrow += 1
						timesum[wrow, :] = Float64[maxtime, NI_now, prop_ever_now, prop_cur_now, NR_now, prop_rec_now, NR_now]
						break
					end

				#	Order infected FIFO by infection_order
					resize!(infected_perm_buf, ninf)
					for i in 1:ninf
						infected_perm_buf[i] = i
					end
					sort!(infected_perm_buf; by = i -> infected_order_buf[i])
					for ii in 1:ninf
						ordered_idx_buf[ii] = infected_idx_buf[infected_perm_buf[ii]]
					end

				#	Process each infected ego
					for p in 1:ninf
						ego_idx = ordered_idx_buf[p]

						#	Neighbors in original column order
							susceptible_neighbor_ids = s_alst_vec[ego_idx]

						#	Ego symptomatic draw (per-day)
							issympt = (rand() < p_symp) ? 1 : 0

						#	Scan neighbors
							for alter_id in susceptible_neighbor_ids
								alter_idx = id_to_idx[alter_id]
								edgwgt    = s_vlst_map[ego_idx][alter_id]

								#	Current cumulative peers-infected
									peersinf = nbrsinf[alter_idx]

								#	Activation logit
									lwact = b_int +
									        (issympt * b_cxn_symp) +
									        (b_close * edgwgt) +
									        (b_cxn_peers * peersinf) +
									        (b_cxn_total * (ninf / (n_nodes / 3))) +
									        (b_cls_x_smp * peersinf * issympt)
									prob_act = exp(lwact) / (1 + exp(lwact))

								#	Activation and transmission
									if rand() < prob_act
										transprob = transmission_method === :weighted ? (1 - (1 - inf_r)^edgwgt) : inf_r
										if rand() < transprob
											#	Update state
												state[alter_idx, INFECTED_COL]         = 1
												state[alter_idx, SUSCEPTIBLE_COL]      = 0
												state[alter_idx, TIME_TO_RECOVERY_COL] = rec_t
												infection_counter += 1
												state[alter_idx, INFECTION_ORDER_COL]  = infection_counter

											#	Remove alter from all susceptible lists
												for k in 1:n_nodes
													remove_id!(s_alst_vec[k], alter_id)
												end

											#	Increment peers-infected for ego’s original neighbors
												for nbr_id in @view alst_vec[ego_idx][2:end]
													nbr_idx = id_to_idx[nbr_id]
													nbrsinf[nbr_idx] += 1
												end
										end
									end
							end

						#	Recovery countdown
							state[ego_idx, TIME_TO_RECOVERY_COL] -= 1
					end

				#	Move recovered
					for i in 1:n_nodes
						if state[i, TIME_TO_RECOVERY_COL] == 0 && state[i, INFECTED_COL] == 1
							state[i, INFECTED_COL]         = 0
							state[i, RECOVERED_COL]        = 1
							state[i, TIME_TO_RECOVERY_COL] = 0
							state[i, INFECTION_ORDER_COL]  = 0
							if immunity_duration !== nothing
								state[i, TIME_TO_RECOVERY_COL] = -immunity_duration
							end
						end
					end

				#	SIRS waning immunity (optional)
					if immunity_duration !== nothing
						for i in 1:n_nodes
							if state[i, TIME_TO_RECOVERY_COL] < 0 && state[i, RECOVERED_COL] == 1
								state[i, TIME_TO_RECOVERY_COL] += 1
							end
						end
						for i in 1:n_nodes
							if state[i, TIME_TO_RECOVERY_COL] == 0 && state[i, RECOVERED_COL] == 1
								state[i, RECOVERED_COL]   = 0
								state[i, SUSCEPTIBLE_COL] = 1
								for j in 1:n_nodes
									if i != j && (unique_ids[i] in alst_vec[j][2:end])
										push!(s_alst_vec[j], unique_ids[i])	# Only in SIRS mode
									end
								end
							end
						end
					end

				#	State consistency check
					for i in 1:n_nodes
						sum_row = state[i, INFECTED_COL] + state[i, SUSCEPTIBLE_COL] + state[i, RECOVERED_COL]
						@assert sum_row == 1 "Invalid state for node $(state[i, ID_COL]): I=$(state[i,2]) S=$(state[i,3]) R=$(state[i,4])"
					end

				#	Record daily metrics
					NI = 0
					NR = 0
					for i in 1:n_nodes
						NI += state[i, INFECTED_COL]
						NR += state[i, RECOVERED_COL]
					end
					prop_ever = (NI + NR) / n_nodes
					prop_cur  = NI / n_nodes
					prop_rec  = (NI + NR) > 0 ? NR / (NI + NR) : 0.0

					wrow += 1
					timesum[wrow, :] = Float64[t, NI, prop_ever, prop_cur, NR, prop_rec, NR]
			end

		#	Stop clock
			total_time = time() - start_time

		#	Trim time series
			inflog = timesum[1:wrow, :]

		#	Return results
			return Dict{String,Any}(
				"infection_log" => inflog,
				"total_time"    => total_time,
				"final_state"   => state,
			)
	end
	@doc raw"""
	**Description**  
	Simulates SIR/SIRS diffusion matching SAS exactly with preserved column order
	for RNG synchronization and cumulative peer effects.

	**Usage**  
	`sirdif(alst, vlst, infectedp, inf_r, rec_t, maxtime, p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp; transmission_method=:weighted, immunity_duration=nothing)`

	**Arguments**
	- `alst`: Adjacency list (Matrix or Vector{Vector})
	- `vlst`: Edge weights aligned to `alst`
	- `infectedp`: Initial infected node IDs
	- `inf_r`: Base transmission probability
	- `rec_t`: Days to recovery (constant)
	- `maxtime`: Simulation horizon (days)
	- `p_symp`: Probability infected is symptomatic
	- `b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp`: Logit coefficients
	- `transmission_method`: `:weighted` or `:simple`
	- `immunity_duration`: Days of immunity after recovery (SIRS mode)

	**Details**
	- Preserves exact column order from adjacency matrix for RNG alignment
	- `nbrsinf` is cumulative and persistent across entire simulation
	- Infected processed in FIFO order by infection time
	- Immediate within-timestep updates when transmission occurs
	- **SAS-style termination**: when infections drop to zero, a single final row is written at `time=maxtime`.

	**Value**
	Dict with `infection_log` (time series), `total_time`, and `final_state`

	**See Also**
	`sim_prep`, `replicate_sas_simulation`
	""" sirdif

#   Exporting Objects
    export arithmetic_mode,
           color_assignment,
           plot_parameter_space,
           sim_parm_space,
           test_sim_parm_space,
		   sim_prep,
		   sirdif

end # module diffustion_sim

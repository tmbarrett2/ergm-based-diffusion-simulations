module diffustion_sim
#   Packages
    using CSV
    using DataFrames
    using RCall
    using Statistics
   
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
			@rput results edgemax b_int b_cxn_symp b_cxn_total b_cxn_peers display
		
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
						X11(width=10, height=12)
				"""
			end
		
		#	Execute R plotting code
			R"""
			# Creating Layout
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
			
			# Looping Over Visualization Matrix
			for (issymptomatic in seq(0, 1, 1)) {
				# Filter for current symptomatic status
				r_panel_data <- results[results$issymptomatic == issymptomatic, ]
				
				for (edgwgt in seq(1, edgemax, 1)) {
					# Isolating Cell Data
					r_cell_data <- r_panel_data[r_panel_data$edgwgt == edgwgt, ]
					
					# Creating Base Plot
					par(mar=c(3, 4, 2, 2))
					plot(NA, type="n", 
					     main = paste('edgwgt =', edgwgt), 
					     xlab = "",
					     ylab = "Probability of Interaction", 
					     family = 'serif',
					     xlim = c(0, 5), 
					     ylim = c(0, 1), 
					     bty = "n", 
					     las = 1)
					
					# Adding Reference Lines
					ref_lines <- seq(0, 1, 0.2)
					for (ref_line in ref_lines) {
						abline(h = ref_line, col = "gray60", lty = 3)
					}
					
					# Plotting X-Axis Label
					mtext(side = 1, text = 'Number of Ill Peers', 
					      col = "black", line = 2, cex = 0.75, family = 'serif')
					
					# Plotting Lines for each ninf level
					ninf_levels <- unique(r_cell_data$ninf)
					for (ninf in ninf_levels) {
						r_line_data <- r_cell_data[r_cell_data$ninf == ninf, ]
						r_line_data <- r_line_data[order(r_line_data$peersinf), ]
						lines(r_line_data$peersinf, r_line_data$prob_act, 
						      lty = 1, col = r_line_data$r_colors[1], lwd = 2)
					}
					
					# Adding Symptomatic label on rightmost panel
					if (edgwgt == edgemax) {
						mtext(side = 4, text = paste("issymptomatic =", issymptomatic), 
						      col = "black", line = 1, cex = 0.75, 
						      family = 'serif', font = 2)
					}
				}
			}
			
			# Adding Plot Title
			par(mar = c(0, 0, 0, 0), bty = 'n')
			plot(0, type = 'n', xlab = ' ', ylab = ' ', 
			     cex.axis = 1.3, xaxt = 'n', yaxt = 'n', 
			     family = 'serif', las = 1, main = ' ')
			text(x = 0.5, y = 0.5, 
			     labels = "SIR Social Feedback Function Plot", 
			     cex = 2.5, family = 'serif', font = 2)
			
			# Adding Color Gradient Legend
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
			
			# Adding Parameter Annotation
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
				paste("b_cls_x_smp =", b_cls_x_smp),
				sep = "  "
			)
			
			text(x = 0.5, y = 0.5, labels = param_text, 
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

#   Exporting Objects
    export arithmetic_mode,
           color_assignment,
           plot_parameter_space,
           sim_parm_space,
           test_sim_parm_space

end # module diffustion_sim

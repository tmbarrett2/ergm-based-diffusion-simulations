#SIR Social Feedback Simulation
#Jonathan H. Morgan, Ph.D., James Moody Ph.D., & Tyler Barret
#27 August 2025

#	To Do:
#	Complete SIRS testing
#	Move Small Speed Test Loop into a Wrapper that R interfaces with
#	Conduct Julia/R Test of Wrapper
#	Conduct Jim's Full Sweep Tests
#	Move Julia Package into an executable
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
using ProgressMeter
using Random
using Statistics
using diffustion_sim

######################
#   TEST FUNCTIONS   #
######################

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
			- inf_r=0.33, rec_t=14, maxtime=200, p_symp=0.5
			- b_int=-0.1, b_close=1.0, b_cxn_peers=-0.5, b_cxn_total=-3.5,
			  b_cxn_symp=-1.5, b_cls_x_smp=-0.1
			- SIRS enabled via immunity_duration=immunity_days
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
				"maxtime"      => 200,		# match R SIR tests
				"p_symp"       => 0.5,
				"b_int"        => -0.1,
				"b_close"      => 1.0,
				"b_cxn_peers"  => -0.5,
				"b_cxn_total"  => -3.5,
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
			final_size = final_row[3]							# cumulative proportion ever infected
			peak_infected = maximum(infection_log[:, 2])
			peak_time = findfirst(x -> x == peak_infected, infection_log[:, 2])

		#	Look for reinfection waves (robust peak detection)
			infected_series = infection_log[:, 2]

			#	Smooth lightly (5-day moving average) to suppress single-day noise
				smoothed = similar(infected_series)
				window = 5
				halfwin = fld(window, 2)
				for i in eachindex(infected_series)
					l = max(1, i - halfwin)
					r = min(length(infected_series), i + halfwin)
					#	mean without importing Statistics
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

			#	Check 1: Conservation (infected + recovered ≤ n)
				total_affected = infection_log[:, 2] .+ infection_log[:, 5]
				checks["conservation"] = all(total_affected .<= n)

			#	Check 2: Reinfection possible after immunity (using filtered peaks)
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
					"duration" => size(infection_log, 1) - 1,
					"n_waves" => n_waves
				),
				"checks" => checks,
				"params" => params,
				"seeds"  => seed_nodes
			)
	end

#	Comparative Test Function
	function compare_sir_sirs(network_data::NamedTuple;
	                         immunity_days::Int=30,
	                         n_runs::Int=10,
	                         verbose::Bool=true)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
			immunity_days::Int: immunity duration for SIRS
			n_runs::Int: number of simulation runs for statistics
			verbose::Bool: print detailed output
		Returns:
			Dict: comparative statistics
		Notes:
			Runs multiple simulations to compare SIR and SIRS dynamics.
		"""
		
		#	Storage for results
			sir_metrics = Dict(
				"final_sizes" => Float64[],
				"peak_infected" => Int[],
				"durations" => Int[]
			)
			
			sirs_metrics = Dict(
				"final_sizes" => Float64[],
				"peak_infected" => Int[],
				"durations" => Int[],
				"n_waves" => Int[]
			)
		
		#	Run simulations
			if verbose
				println("="^60)
				println("Running Comparative Tests")
				println("="^60)
				println("Number of runs: $n_runs")
			end
			
			for run in 1:n_runs
				#	SIR simulation
					sir_result = test_sir_model(network_data; 
					                           seed=run * 100, 
					                           verbose=false)
					push!(sir_metrics["final_sizes"], sir_result["metrics"]["final_size"])
					push!(sir_metrics["peak_infected"], sir_result["metrics"]["peak_infected"])
					push!(sir_metrics["durations"], sir_result["metrics"]["duration"])
				
				#	SIRS simulation
					sirs_result = test_sirs_model(network_data;
					                             immunity_days=immunity_days,
					                             seed=run * 100,
					                             verbose=false)
					push!(sirs_metrics["final_sizes"], sirs_result["metrics"]["final_size"])
					push!(sirs_metrics["peak_infected"], sirs_result["metrics"]["peak_infected"])
					push!(sirs_metrics["durations"], sirs_result["metrics"]["duration"])
					push!(sirs_metrics["n_waves"], sirs_result["metrics"]["n_waves"])
					
				if verbose && run % 5 == 0
					println("  Completed run $run/$n_runs")
				end
			end
		
		#	Calculate statistics
			if verbose
				println("\n" * "="^40)
				println("Comparative Results ($n_runs runs):")
				println("\nSIR Model:")
				println("  Final size: $(round(mean(sir_metrics["final_sizes"]) * 100, digits=1))% " *
				       "(SD: $(round(std(sir_metrics["final_sizes"]) * 100, digits=1))%)")
				println("  Peak infected: $(round(mean(sir_metrics["peak_infected"]), digits=1)) " *
				       "(SD: $(round(std(sir_metrics["peak_infected"]), digits=1)))")
				println("  Duration: $(round(mean(sir_metrics["durations"]), digits=1)) days " *
				       "(SD: $(round(std(sir_metrics["durations"]), digits=1)))")
				
				println("\nSIRS Model (immunity = $immunity_days days):")
				println("  Final size: $(round(mean(sirs_metrics["final_sizes"]) * 100, digits=1))% " *
				       "(SD: $(round(std(sirs_metrics["final_sizes"]) * 100, digits=1))%)")
				println("  Peak infected: $(round(mean(sirs_metrics["peak_infected"]), digits=1)) " *
				       "(SD: $(round(std(sirs_metrics["peak_infected"]), digits=1)))")
				println("  Duration: $(round(mean(sirs_metrics["durations"]), digits=1)) days " *
				       "(SD: $(round(std(sirs_metrics["durations"]), digits=1)))")
				println("  Average waves: $(round(mean(sirs_metrics["n_waves"]), digits=1))")
			end
		
		#	Return comparative results
			return Dict(
				"sir" => sir_metrics,
				"sirs" => sirs_metrics
			)
	end

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
    sirs_results = test_sirs_model(network_data; immunity_days=30, seed=42, verbose=true)

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


#   Comparisative Tests



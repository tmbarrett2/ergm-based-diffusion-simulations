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
using RCall
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

#	Helper Function for SIR Comparator
	function summarize_ts(df::DataFrame)
		peak_val = maximum(df.Number_Infected)
		peak_idx = findfirst(==(peak_val), df.Number_Infected)
		Dict(
			:final_size_prop => df.ProportionEverInfected[end],
			:peak_infected   => peak_val,
			:peak_time       => df.Time[peak_idx],
			:duration_days   => df.Time[end] - df.Time[1]
		)
	end

#	Helper: resolve a column by case-insensitive candidate names
	function resolve_col(df::DataFrame, candidates::Vector{Symbol})
		"""Resolve a DataFrame column (case-insensitive) from a list of candidate names."""
		name_map = Dict{String,Symbol}()
		for nm in names(df)
			name_map[strip(lowercase(String(nm)))] = nm
		end
		for cand in candidates
			key = strip(lowercase(String(cand)))
			if haskey(name_map, key)
				return name_map[key]
			end
		end
		return nothing
	end

#	Helper: compute basic summary metrics for a normalized frame
	function compute_summary_metrics(df::DataFrame)
		"""Return Dict with :final_size, :peak_infected, :peak_time, :duration for the input series."""
		peak_infected = maximum(df.Number_Infected)
		peak_time_idx = findfirst(x -> x == peak_infected, df.Number_Infected)
		return Dict(
			:final_size => df.ProportionEverInfected[end],
			:peak_infected => peak_infected,
			:peak_time => df.Time[peak_time_idx],
			:duration => df.Time[end] - df.Time[1]
		)
	end

#	Comparator: Julia SIR vs SAS CSV (robust CSV column resolution)
	function compare_sir_sim_to_sas(sim_or_test::Dict, sas_csv::String;
	                                series::Symbol=:NI,
	                                sas_case::Union{Nothing,Int}=nothing,
	                                show_plot::Bool=true,
	                                save_path::Union{Nothing,String}=nothing,
	                                display::String=":0")
		"""Args:
			sim_or_test::Dict: either sirdif() result or test_*_model() result
			sas_csv::String: path to SAS CSV file
			series::Symbol: :NI, :prop_ever, :prop_cur, :NR, or :prop_rec
			sas_case::Union{Nothing,Int}: specific case number to extract from SAS data
			show_plot::Bool: whether to display plot interactively
			save_path::Union{Nothing,String}: optional PDF output path
			display::String: X11 display for plotting
		Returns:
			Dict: comparison metrics and aligned data
		Notes:
			Compares Julia simulation with SAS output.
			Handles both raw sirdif output and test wrapper results.
			Infers network size from initial conditions.
		"""
		
		#	Unwrap test wrapper if present
			sim_results = haskey(sim_or_test, "results") ? sim_or_test["results"] : sim_or_test
			
		#	Validate infection log exists
			if !haskey(sim_results, "infection_log")
				error("Input must contain 'infection_log' key")
			end
		
		#	Extract Julia infection log
			log = sim_results["infection_log"]
			
		#	Validate log structure
			if size(log, 2) < 6
				error("infection_log must have at least 6 columns")
			end
			
		#	Convert to DataFrame
			df_julia = DataFrame(
				Time = Int.(log[:, 1]),
				Number_Infected = Int.(log[:, 2]),
				ProportionEverInfected = log[:, 3],
				ProportionCurrentlyInfected = log[:, 4],
				Number_Recovered = Int.(log[:, 5]),
				ProportionRecovered = log[:, 6]
			)
		
		#	Infer network size from initial conditions
			n0_infected = df_julia.Number_Infected[1]
			n0_prop_ever = df_julia.ProportionEverInfected[1]
			if n0_prop_ever <= 0
				error("Cannot infer network size: initial proportion ever infected is zero")
			end
			n_estimated = round(Int, n0_infected / n0_prop_ever)
			if n_estimated <= 0
				error("Inferred network size must be positive")
			end
		
		#	Load SAS CSV
			df_sas_raw = CSV.read(sas_csv, DataFrame)

		#	Resolve required/optional columns (case-insensitive)
			col_TIME = resolve_col(df_sas_raw, Symbol.(["TIME","Time","time"]))
			col_EVER = resolve_col(df_sas_raw, Symbol.(["PROPEVERINF","PropEverInf","propeverinf"]))
			col_CURR = resolve_col(df_sas_raw, Symbol.(["PROPCURRINF","PropCurrInf","propcurrinf"]))
			col_REC  = resolve_col(df_sas_raw, Symbol.(["PROPREC","PropRec","proprec"]))
			col_CASE = resolve_col(df_sas_raw, Symbol.(["case","CASE","Case"]))

			if col_TIME === nothing
				error("SAS CSV missing TIME column. Found columns: $(names(df_sas_raw))")
			end
			if col_EVER === nothing
				error("SAS CSV missing PROPEVERINF column. Found columns: $(names(df_sas_raw))")
			end

		#	Optional filter by case
			if sas_case !== nothing && col_CASE !== nothing
				df_sas_raw = filter(row -> row[col_CASE] == sas_case, df_sas_raw)
			elseif col_CASE !== nothing
				ucases = unique(df_sas_raw[!, col_CASE])
				if length(ucases) > 1
					df_sas_raw = filter(row -> row[col_CASE] == first(ucases), df_sas_raw)
				end
			end

		#	Standardize schema
			rename!(df_sas_raw, Dict(col_TIME => :Time, col_EVER => :ProportionEverInfected))
			if col_CURR !== nothing
				rename!(df_sas_raw, Dict(col_CURR => :ProportionCurrentlyInfected))
			else
				df_sas_raw.ProportionCurrentlyInfected = fill(NaN, nrow(df_sas_raw))
			end

			if col_REC !== nothing
				rename!(df_sas_raw, Dict(col_REC => :ProportionRecovered))
			else
				df_sas_raw.ProportionRecovered = fill(NaN, nrow(df_sas_raw))
			end

		#	Back-fill missing columns from identities
			E  = df_sas_raw.ProportionEverInfected
			Rf = df_sas_raw.ProportionRecovered
			Pc = df_sas_raw.ProportionCurrentlyInfected
			for i in 1:nrow(df_sas_raw)
				if (isnan(Rf[i]) || ismissing(Rf[i])) && !(isnan(Pc[i]) || ismissing(Pc[i])) && E[i] > 0
					Rf[i] = 1 - (Pc[i] / E[i])
				end
			end
			Rf .= coalesce.(Rf, 0.0)
			for i in 1:nrow(df_sas_raw)
				if isnan(Pc[i]) || ismissing(Pc[i])
					Pc[i] = E[i] * (1 - Rf[i])
				end
			end

		#	Finalize SAS frame
			df_sas = DataFrame(
				Time = Int.(df_sas_raw.Time),
				ProportionEverInfected = E,
				ProportionCurrentlyInfected = Pc,
				ProportionRecovered = Rf
			)
			df_sas.Number_Infected = round.(Int, n_estimated .* df_sas.ProportionCurrentlyInfected)
			df_sas.Number_Recovered = round.(Int, n_estimated .* df_sas.ProportionEverInfected .* df_sas.ProportionRecovered)
		
		#	Define series mapping
			series_map = Dict(
				:NI => (:Number_Infected, "Number Infected"),
				:prop_ever => (:ProportionEverInfected, "Proportion Ever Infected"),
				:prop_cur => (:ProportionCurrentlyInfected, "Proportion Currently Infected"),
				:NR => (:Number_Recovered, "Number Recovered"),
				:prop_rec => (:ProportionRecovered, "Proportion Recovered")
			)
			if !haskey(series_map, series)
				error("Invalid series: $series. Choose from: $(keys(series_map))")
			end
			series_col, series_label = series_map[series]
		
		#	Align time series
			julia_series = select(df_julia, [:Time, series_col])
			rename!(julia_series, series_col => :julia)
			sas_series = select(df_sas, [:Time, series_col])
			rename!(sas_series, series_col => :sas)
			aligned = innerjoin(julia_series, sas_series, on=:Time)
			sort!(aligned, :Time)
			n_pairs = nrow(aligned)
			if n_pairs == 0
				error("No overlapping time points between Julia and SAS")
			end
		
		#	Error metrics
			correlation = n_pairs > 1 ? cor(aligned.julia, aligned.sas) : NaN
			mae = mean(abs.(aligned.julia .- aligned.sas))
			rmse = sqrt(mean((aligned.julia .- aligned.sas).^2))
			errors = Dict(
				:n => n_pairs,
				:correlation => correlation,
				:mae => mae,
				:rmse => rmse
			)
		
		#	Summary metrics
			julia_metrics = compute_summary_metrics(df_julia)
			sas_metrics = compute_summary_metrics(df_sas)
		
		#	PLOT (RCall)
			@rput aligned series_label show_plot save_path display correlation mae rmse
			R"""
			if (!is.null(save_path)) {
				pdf(save_path, width=8, height=5)
				plot_device <- "pdf"
			} else if (show_plot) {
				Sys.setenv(DISPLAY = display)
				X11(type="cairo", width=8, height=5)
				plot_device <- "x11"
			} else {
				pdf(tempfile(fileext=".pdf"), width=8, height=5)
				plot_device <- "temp"
			}
			par(mar=c(4,4,2,1))
			plot(aligned$Time, aligned$sas,
			     type="l", lwd=2, col="gray40",
			     xlab="Time (days)", ylab=series_label,
			     las=1, bty="n",
			     main="SIR Model Comparison: Julia vs SAS")
			lines(aligned$Time, aligned$julia, lwd=2, col="black")
			legend("topleft", legend=c("SAS","Julia"),
			       col=c("gray40","black"), lwd=2, bty="n")
			mtext(paste("Correlation:", round(correlation, 3),
			            " MAE:", round(mae, 2),
			            " RMSE:", round(rmse, 2)),
			      side=3, line=0.5, cex=0.8)
			if (plot_device %in% c("pdf","temp")) dev.off()
			"""
			plot_saved = save_path !== nothing
		
		#	Return comparison results
			return Dict(
				"julia_metrics" => julia_metrics,
				"sas_metrics" => sas_metrics,
				"errors" => errors,
				"aligned" => aligned,
				"plot_saved" => plot_saved,
				"network_size" => n_estimated
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
    sirs_results = test_sirs_model(network_data; immunity_days=14, seed=42, verbose=true)

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
	comp_results = compare_sir_sim_to_sas(sir_results["results"],
										  "/workspace/data/sim_test_data/SAS_Simulation_Outputs.csv";
										  series     = :NI,                 # choose :NI, :prop_ever, :prop_cur, :NR
										  show_plot  = true,                # open X11 window
										  save_path  = nothing,             # or provide a PDF path like "sir_overlay.pdf"
										  display    = ":100"               # adjust if remote DISPLAY differs
					)

#	COME BACK HERE!!!!!
#	ERROR: MethodError: Cannot `convert` an object of type String to an object of type Symbol
# 	The function `convert` exists, but no method is defined for this combination of argument types.
#	Should use parse() here.

# 	Inspect results
	println("Comparison errors: ", comp_results[:errors])
	println("Julia metrics: ", comp_results[:julia_metrics])
	println("SAS metrics: ", comp_results[:sas_metrics])


#SIR Social Feedback Simulation
#Jonathan H. Morgan, Ph.D., James Moody Ph.D., & Tyler Barret
#27 August 2025

#	To Do:
#	Complete Simulation Test: Look at specific run to see why SAS is consistly scoring a higher infection rate.
#	Move Small Speed Test Loop into a Wrapper that R interfaces with
#	Conduct Julia/R Test of Wrapper
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
				:b_int        => -0.1,
				:b_close      =>  1.0,
				:b_cxn_peers_root  => -0.5,
				:b_cxn_global_root => -3.5,
				:b_cxn_symp_root   => -1.5,
				:b_cls_x_smp       => -0.1,
				:p_symp            =>  0.75,
				:progress_desc     => "SAS IML-style simulation (design-driven)"
			)

		#	Build coefficient design (4×4×6 = 96 rows)
			b_peer_v, b_glob_v, b_self_v = replicate_sas_simulation_build_design(
				params[:b_cxn_peers_root], params[:b_cxn_global_root], params[:b_cxn_symp_root]
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

		# 	coefficients (varied) + constants (for traceability)
			b_cxn_peer_v   = Vector{Float64}(undef, max_out)
			b_cxn_global_v = Vector{Float64}(undef, max_out)
			b_self_v_out   = Vector{Float64}(undef, max_out)
			b_int_v        = Vector{Float64}(undef, max_out)
			b_close_v      = Vector{Float64}(undef, max_out)
			b_cls_x_smp_v  = Vector{Float64}(undef, max_out)
			p_symp_v       = Vector{Float64}(undef, max_out)
			beta_v         = Vector{Float64}(undef, max_out)

		# 	outcomes (per time step)
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
					# Random single seed (match SAS sampling of one node)
					seednode = rand(1:n_nodes)
					infectedp = [seednode]

					# Run sirdif with this row's coefficients
					result = sirdif(
						alst, vlst, infectedp,
						inf_r, rec_t, maxtime,
						params[:p_symp],
						params[:b_int], params[:b_close],
						b_peers, b_glob, b_self, params[:b_cls_x_smp];
						transmission_method = :weighted,
						immunity_duration   = nothing
					)

					# Copy trajectory rows
					logmat = result["infection_log"]::Matrix{Float64}
					T      = size(logmat, 1)

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

						time_v[w]         = logmat[k, 1]
						propeverinf_v[w]  = logmat[k, 3]
						propcurrinf_v[w]  = logmat[k, 4]
						proprec_v[w]      = logmat[k, 6]
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
				proprec      = proprec_v[1:last]
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
		outcome::Symbol = :propeverinf,                # :propeverinf | :propcurrinf | :proprec
		display::String = ":100",
		save_path::Union{String,Nothing} = nothing,
		show_plot::Bool = true,
		seed::Union{Int,Nothing} = nothing
	)
		"""
		Args:
			alst, vlst: Network adjacency and weights (matrix or vector-of-vectors; same form).
			inf_r::Float64: Base transmission probability per contact (β).
			rec_t::Int: Recovery time in days.
			maxtime::Int: Maximum simulation days (≥ 1).
			num_iter::Int: Iterations per coefficient setting (≥ 1).
			sas_csv_path::AbstractString: Path to SAS PROC IML output CSV.
		Returns:
			DataFrame: One per-time median comparison table with constants → coefficients → SAS medians → Julia medians → absolute deltas.
		Notes:
			- Validates exact match of (time, b_cxn_peer, b_cxn_global, b_self) grids.
			- Facets in R: rows = Global effect, cols = Peer effect, line color = Symptomatic.
			- Plot uses Julia solid lines vs SAS dashed lines.
			- Bottom annotation shows constants: β, p_symp, b_int, b_close, b_cls_x_smp, rec_t, maxtime.
		"""

		#	Seed (optional)
			if seed !== nothing
				Random.seed!(seed)
			end

		#	Read SAS CSV + aggregates
			df_sas, med_time_sas, _final_per_run_sas, _med_final_sas = sas_simulation_aggregator(sas_csv_path)

		#	Matching Parameters
			num_iter_sas  = length(unique(df_sas.itter))
			maxtime_sas   = maximum(df_sas.time)
		
		#	Run Julia simulation (shows its own progress bar)
			df_julia = replicate_sas_simulation(alst, vlst, inf_r, rec_t, Int(maxtime_sas), num_iter_sas)

		#	Aggregate Julia medians/finals 
			julia_summ = simulation_comparer_build_julia_medians(df_julia)
			med_time_julia = julia_summ.med_time	# :time, b_* + *_med_julia

		#	Harmonize SAS median column names → *_med_sas
			rename!(med_time_sas, Dict(
				:propeverinf_med => :propeverinf_med_sas,
				:propcurrinf_med => :propcurrinf_med_sas,
				:proprec_med     => :proprec_med_sas
			))

		#	Validate coefficients only (not time grid)
			coef_keys = [:b_cxn_peer, :b_cxn_global, :b_self]
			levels_sas   = sort(unique(eachrow(med_time_sas[:, coef_keys])))
			levels_julia = sort(unique(eachrow(med_time_julia[:, coef_keys])))
			@assert levels_sas == levels_julia "Coefficient sets differ."

		#	Join on (time + coefficients) using intersection of times
			comp = innerjoin(
				med_time_sas, med_time_julia,
				on = [:time, :b_cxn_peer, :b_cxn_global, :b_self],
				makeunique = true
			)

		#	Constants used in Julia run (held fixed)
			const_p_symp      = 0.75
			const_b_int       = -0.1
			const_b_close     =  1.0
			const_b_cls_x_smp = -0.1

		#	Absolute deltas: abs(SAS − Julia)
			propeverinf_delta = abs.(comp.propeverinf_med_sas .- comp.propeverinf_med_julia)
			propcurrinf_delta = abs.(comp.propcurrinf_med_sas .- comp.propcurrinf_med_julia)
			proprec_delta     = abs.(comp.proprec_med_sas     .- comp.proprec_med_julia)

		#	Assemble comparison table (pre-allocated constant columns first)
			nrows       = nrow(comp)
			beta_v      = fill(inf_r, nrows)
			p_symp_v    = fill(const_p_symp, nrows)
			b_int_v     = fill(const_b_int, nrows)
			b_close_v   = fill(const_b_close, nrows)
			b_cls_v     = fill(const_b_cls_x_smp, nrows)
			rec_t_v     = fill(rec_t, nrows)
			maxtime_v   = fill(maxtime, nrows)

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
				proprec_delta     = proprec_delta
			)

		#	Plot in R: rows = Global effect, cols = Peer effect, color = Symptomatic
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

		#	Open graphics device if saving
			if save_path !== nothing
				R"""
				pdf($save_path, width=10, height=12)
				"""
			elseif show_plot
				R"""
				#	Set DISPLAY
					Sys.setenv(DISPLAY=':100')
					Sys.getenv("DISPLAY")      
					
				#	Creating X11 Window
					X11(type = "cairo", width=10, height=12)
				"""
			end

		#	Plotting
			R"""
			# ---- Julia->R injected args as local vars ----
			col_j <- $ycol_j
			col_s <- $ycol_s
			ylab  <- $y_label
			ttl   <- $title_label
			subt  <- $subtitle_label
			param <- $param_text
			spath <- $save_path

			# ---- Data prep ----
			df <- comp_table
			df[["b_cxn_peer"]]   <- as.numeric(df[["b_cxn_peer"]])
			df[["b_cxn_global"]] <- as.numeric(df[["b_cxn_global"]])
			df[["b_self"]]       <- as.numeric(df[["b_self"]])
			df[["time"]]         <- as.numeric(df[["time"]])

			# 	Facet levels
			peer_levels <- sort(unique(df[["b_cxn_peer"]]))
			glob_levels <- sort(unique(df[["b_cxn_global"]]))
			self_levels <- sort(unique(df[["b_self"]]))

			# 	Layout: rows = Global effect, cols = Peer effect
			layout(matrix(1:(length(glob_levels)*length(peer_levels)),
							nrow=length(glob_levels), byrow=TRUE))
			par(family="serif", mar=c(3.5,4,2.5,1.5), las=1)

			# 	Ticks helper (Dataplot-style)
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

			# 	Color palette for Symptomatic lines
			pal_base <- c("#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#4477AA", "#BBBBBB")
			if (length(self_levels) <= length(pal_base)) {
			pal <- pal_base[seq_along(self_levels)]
			} else {
			pal <- grDevices::colorRampPalette(pal_base)(length(self_levels))
			}

			# 	Panels
				for (gi in seq_along(glob_levels)) {
					gval <- glob_levels[gi]
					for (pi in seq_along(peer_levels)) {
						# 	Slice panel data
							pval <- peer_levels[pi]
							sub  <- df[df[["b_cxn_global"]] == gval & df[["b_cxn_peer"]] == pval, ]

							ylim <- range(c(sub[[ col_j ]], sub[[ col_s ]]), na.rm=TRUE)
							ylim[1] <- min(0, ylim[1])

						#	Base plot
							plot(NA, type="n",
								xlim=range(sub[["time"]], na.rm=TRUE), ylim=ylim,
								xlab="Time", ylab=ylab, tck=0.015, xaxt='n', bty='L', las=1,
								main=paste0("Peer effect = ", sprintf("%.2f", pval),
											"   |   Global effect = ", sprintf("%.2f", gval)))
							axis(1, padj=0.75, tck=0.015)
							dataplot_tick_function()

						# 	Lines by symptomatic (self) level
							for (si in seq_along(self_levels)) {
								#	Isolating Elements
									sval <- self_levels[si]
									ss <- sub[sub[["b_self"]] == sval, ]
									ss <- ss[order(ss[["time"]]), ]

								# 	Julia (solid)
									lines(ss[["time"]], ss[[ col_j ]], col=pal[si], lwd=2, lty=1)

								# 	SAS (dashed)
									lines(ss[["time"]], ss[[ col_s ]], col=pal[si], lwd=2, lty=2)
							}

							# Legend on last panel
							#if (gi == length(glob_levels) && pi == length(peer_levels)) {
							#leg_lbls <- c(paste0("Symptomatic = ", sprintf("%.2f", self_levels)),
							#				"Julia median", "SAS median")
							#leg_col  <- c(pal, "black", "black")
							#leg_lty  <- c(rep(1, length(self_levels)), 1, 2)
							#leg_lwd  <- c(rep(2, length(self_levels)), 2, 2)
							#legend("topright", legend=leg_lbls, col=leg_col, lty=leg_lty, lwd=leg_lwd,
							#		cex=0.85, bty="n")
							#}
					}
				}

			# 	Title row
				#par(mar = c(0, 0, 0, 0), bty = "n")
				#plot(0, type = "n", xlab = " ", ylab = " ", xaxt = "n", yaxt = "n", main = " ")
				#text(x = 0.5, y = 0.7, labels = ttl, cex = 1.8, family = "serif", font = 2)
				#text(x = 0.5, y = 0.4, labels = subt, cex = 1.0, family = "serif")

				# Bottom annotation of constants
				#par(mar = c(0, 0, 0, 0), bty = "n")
				#plot(0, type = "n", xlab = " ", ylab = " ", xaxt = "n", yaxt = "n", main = " ")
				#text(x = 0.5, y = 0.5, labels = param, cex = 1.1, family = "serif")

			#	Saving
			if (!is.null(spath)) dev.off()
			"""

		#	Assembling Result
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

#   Comparisative Tests: Replication Analysis

#   Load network from GraphML file
#	network_data = sim_prep("/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml", sas_transformation=false)
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

#	Build 5-Node Line Network
	function build_line5()
		"""
		Args:
			None
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (alst, vlst)
		Notes:
			Matrix format with id in col 1; neighbor ids in cols 2..4 (0 = none).
			Weights set to 1.0 for each existing edge; zeros elsewhere.
		"""
		#	Allocate adjacency with id column
			alst = zeros(Int, 5, 4)
			alst[:, 1] = 1:5

		#	Wire neighbors (1–2–3–4–5)
			alst[1, 2:4] = [2, 0, 0]
			alst[2, 2:4] = [1, 3, 0]
			alst[3, 2:4] = [2, 4, 0]
			alst[4, 2:4] = [3, 5, 0]
			alst[5, 2:4] = [4, 0, 0]

		#	Weights: 1.0 where an edge exists (neighbor cols only)
			vlst = zeros(Float64, size(alst))
			vlst[:, 2:4] .= Float64.(alst[:, 2:4] .> 0)

		#	Return
			return alst, vlst
	end

#	Build 5-Node Star Network (center=1)
	function build_star5()
		"""
		Args:
			None
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (alst, vlst)
		Notes:
			Undirected star encoded row-wise; center 1 connected to {2,3,4,5}.
		"""
		#	Allocate adjacency with id column
			alst = zeros(Int, 5, 5)          # center has 4 alters
			alst[:, 1] = 1:5

		#	Wire neighbors
			alst[1, 2:5] = [2, 3, 4, 5]      # vector, not tuple
			alst[2, 2]   = 1
			alst[3, 2]   = 1
			alst[4, 2]   = 1
			alst[5, 2]   = 1

		#	Weights aligned to alst (1.0 on edges in neighbor cols)
			vlst = zeros(Float64, size(alst))
			vlst[:, 1] = Float64.(alst[:, 1])        # carry id in col 1 (for consistency)
			vlst[:, 2:end] .= Float64.(alst[:, 2:end] .> 0)

		#	Return
			return alst, vlst
	end

#	Tier 0 — Test 1: All-On Per-Contact
	function test_t0_all_on(; seed::Int=1234, maxtime::Int=10, rec_t::Int=14)
		"""
		Args:
			seed::Int: RNG seed (for completeness)
			maxtime::Int: simulation horizon (days)
			rec_t::Int: days to recovery
		Returns:
			NamedTuple: (
				line_log, star_log, line_out, star_out, params
			)
		Notes:
			Forces activation ~ 1 via large intercept and sets per-contact transmission = 1
			(`transmission_method = :simple`). Should yield near-deterministic spread governed
			primarily by topology.
		"""
		#	Seed RNG
			Random.seed!(seed)
		#	Build networks
			line_alst, line_vlst = build_line5()
			star_alst, star_vlst = build_star5()
		#	Params: activation on; transmission always succeeds
			p_symp       = 0.0
			b_int        = 20.0         # logistic(20) ≈ 1
			b_close      = 0.0
			b_cxn_peers  = 0.0
			b_cxn_total  = 0.0
			b_cxn_symp   = 0.0
			b_cls_x_smp  = 0.0
			inf_r        = 1.0
			tx_method    = :simple
		#	Run: line (seed at node 1)
			line_out = sirdif(
				line_alst, line_vlst,
				[1], inf_r, rec_t, maxtime,
				p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp;
				transmission_method = tx_method,
				immunity_duration   = nothing
			)
			line_log = line_out["infection_log"]
		#	Run: star (seed at node 1, the center)
			star_out = sirdif(
				star_alst, star_vlst,
				[1], inf_r, rec_t, maxtime,
				p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp;
				transmission_method = tx_method,
				immunity_duration   = nothing
			)
			star_log = star_out["infection_log"]
		#	Light invariants
			@assert size(line_log, 2) == 7 "line_log must have 7 columns"
			@assert size(star_log, 2) == 7 "star_log must have 7 columns"
			@assert line_log[1, 1] == 0.0 "first row is time 0 for line"
			@assert star_log[1, 1] == 0.0 "first row is time 0 for star"
		#	Return artifacts and params
			return (
				line_log = line_log,
				star_log = star_log,
				line_out = line_out,
				star_out = star_out,
				params = (
					p_symp      = p_symp,
					b_int       = b_int,
					b_close     = b_close,
					b_cxn_peers = b_cxn_peers,
					b_cxn_total = b_cxn_total,
					b_cxn_symp  = b_cxn_symp,
					b_cls_x_smp = b_cls_x_smp,
					inf_r       = inf_r,
					tx_method   = tx_method,
					rec_t       = rec_t,
					maxtime     = maxtime,
					seed        = seed
				)
			)
	end

#	Test 1
	res = test_t0_all_on(seed=1234, maxtime=10, rec_t=14)

#	Quick peek
	first(res.line_log, 6), first(res.star_log, 6)

#	Access pieces
	line_log = res.line_log
	star_log = res.star_log
	params   = res.params



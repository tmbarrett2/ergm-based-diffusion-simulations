#Simulator Tests: UT-1.6 and UT-1.7
#Jonathan H. Morgan, Ph.D.
#29 September 2026

#   Activating Local Environment
	using Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	Pkg.status()

################
#   PACKAGES   #
################
using DataFrames
using Random
using Statistics
using diffusion_sim

#   Frozen Copy of the Validated Simulator
	include(joinpath(@__DIR__, "reference_sirdif.jl"))

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Simulator Tests: build a small network in sim_prep format
	function build_test_network(neighbors::Vector{Vector{Int}})
		"""
		Args:
			neighbors::Vector{Vector{Int}}: neighbors[i] lists the IDs tied to node i (IDs are 1:n)
		Returns:
			NamedTuple: (alst::Matrix{Int}, vlst::Matrix{Float64}) in sim_prep format
		Notes:
			Column 1 holds the node ID; neighbor slots are padded with 0.
			All edge weights are 1.0. Ties must be listed in both directions.
		"""

		#	Size the matrices
			n = length(neighbors)
			max_degree = max(1, maximum(length.(neighbors)))
			alst = zeros(Int, n, max_degree + 1)
			vlst = zeros(Float64, n, max_degree + 1)

		#	Fill IDs, neighbors, and weights
			for i in 1:n
				alst[i, 1] = i
				vlst[i, 1] = Float64(i)
				for (col, j) in enumerate(neighbors[i])
					alst[i, col + 1] = j
					vlst[i, col + 1] = 1.0
				end
			end

		#	Assembling Result
			return (alst = alst, vlst = vlst)
	end

#	Helper Function for test_waning_schedule: expected daily counts on the path network
	function expected_path_counts(infection_days::Vector{Int}, D::Int, M::Int, t_end::Int)
		"""
		Args:
			infection_days::Vector{Int}: day each node is infected (row index in the log)
			D::Int: fixed recovery duration in days
			M::Int: fixed immunity duration in days
			t_end::Int: last day to tabulate
		Returns:
			DataFrame: time, n_infected, n_recovered, n_ever for days 0:t_end
		Notes:
			A node infected on day e is infected in rows e to e+D-1, recovers at
			the end of day e+D, is recovered in rows e+D to e+D+M-1, and becomes
			susceptible at the end of day e+D+M.
		"""

		#	Allocate counts
			days = collect(0:t_end)
			n_inf = zeros(Int, length(days))
			n_rec = zeros(Int, length(days))
			n_ever = zeros(Int, length(days))

		#	Tabulate each node's schedule
			for e in infection_days
				r = e + D
				for (k, t) in enumerate(days)
					if e <= t < r
						n_inf[k] += 1
					end
					if r <= t < r + M
						n_rec[k] += 1
					end
					if t >= e
						n_ever[k] += 1
					end
				end
			end

		#	Assembling Result
			return DataFrame(time = days, n_infected = n_inf, n_recovered = n_rec, n_ever = n_ever)
	end

#	Duration Draw Properties
	function test_draw_duration(; n_draws::Int = 100_000, seed::Int = 2026, verbose::Bool = true)
		"""
		Args:
			n_draws::Int: number of Weibull draws per check (default = 100_000)
			seed::Int: seed for the explicit generator (default = 2026)
			verbose::Bool: print a report (default = true)
		Returns:
			NamedTuple: (checks::Dict{String,Bool}, all_passed::Bool)
		Notes:
			Checks the one-day floor, the Weibull mean, that FixedDuration makes
			no random draws, and that invalid specifications are rejected.
		"""

		#	Initialize checks
			checks = Dict{String, Bool}()

		#	Floor: a scale below one day still yields at least one day
			rng = Xoshiro(seed)
			short = [draw_duration(WeibullDuration(0.5, 2.0), rng) for _ in 1:n_draws]
			checks["floor_one_day"] = minimum(short) >= 1

		#	Mean: scale 10, shape 2 has mean 10 * Γ(1.5) = 5 * sqrt(pi)
			rng = Xoshiro(seed + 1)
			long = [draw_duration(WeibullDuration(10.0, 2.0), rng) for _ in 1:n_draws]
			target = 5.0 * sqrt(pi)
			checks["weibull_mean"] = abs(mean(long) - target) < 0.10

		#	FixedDuration leaves the random number stream untouched
			rng_a = Xoshiro(seed + 2)
			rng_b = copy(rng_a)
			fixed_days = draw_duration(FixedDuration(14), rng_a)
			checks["fixed_value"] = fixed_days == 14
			checks["fixed_no_draw"] = rand(rng_a) == rand(rng_b)

		#	Invalid specifications are rejected
			checks["reject_fixed_zero"] = try
				FixedDuration(0)
				false
			catch e
				e isa DomainError
			end
			checks["reject_weibull_scale"] = try
				WeibullDuration(-1.0, 2.0)
				false
			catch e
				e isa DomainError
			end
			checks["reject_weibull_shape"] = try
				WeibullDuration(10.0, 0.0)
				false
			catch e
				e isa DomainError
			end

		#	Report
			all_passed = all(values(checks))
			if verbose
				println("=" ^ 60)
				println("Duration Draw Properties")
				println("=" ^ 60)
				println("  Minimum draw at scale 0.5: $(minimum(short))")
				println("  Mean at scale 10, shape 2: $(round(mean(long), digits=4))  (target $(round(target, digits=4)))")
				for (name, passed) in sort(collect(checks))
					println("  $(rpad(name, 24)) $(passed ? "PASS" : "FAIL")")
				end
				println("RESULT: $(all_passed ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")")
			end

		#	Assembling Result
			return (checks = checks, all_passed = all_passed)
	end
	@doc raw"""
	**Description**
	Check the properties of `draw_duration`: the one-day floor, the Weibull mean, that `FixedDuration` consumes no random numbers, and that invalid specifications throw `DomainError`.

	**Usage**
	`test_draw_duration(; n_draws=100_000, seed=2026, verbose=true)`

	**Value**
	A `NamedTuple` with `checks::Dict{String,Bool}` and `all_passed::Bool`.

	**See Also**
	`draw_duration`, `FixedDuration`, `WeibullDuration`
	""" test_draw_duration

#	Weibull Recovery Timing on an Isolated Node
	function test_weibull_recovery_timing(; scales::Vector{Float64} = [0.5, 10.0],
											shapes::Vector{Float64} = [0.8, 2.0],
											seeds = 1:50, maxtime::Int = 400,
											verbose::Bool = true)
		"""
		Args:
			scales::Vector{Float64}: Weibull scales to test (default = [0.5, 10.0])
			shapes::Vector{Float64}: Weibull shapes to test (default = [0.8, 2.0])
			seeds: seeds for the explicit generator (default = 1:50)
			maxtime::Int: simulation horizon (default = 400)
			verbose::Bool: print a report (default = true)
		Returns:
			NamedTuple: (results::DataFrame, all_passed::Bool)
		Notes:
			A single isolated seed makes the recovery draw the first use of rng,
			so the expected duration is draw_duration(spec, copy(rng)). The
			node must be infected in rows 0 to D-1 and recovered in row D.
		"""

		#	Build a one-node network
			net = build_test_network([Int[]])

		#	Pre-allocate results
			n_cases = length(scales) * length(shapes) * length(seeds)
			results = DataFrame(scale = zeros(n_cases), shape = zeros(n_cases), seed = zeros(Int, n_cases),
								drawn_days = zeros(Int, n_cases), passed = falses(n_cases))

		#	Iterate over specifications and seeds
			k = 0
			for scale in scales, shape in shapes, s in seeds
				#	Predict the recovery day from a copy of the generator
					k += 1
					spec = WeibullDuration(scale, shape)
					D = draw_duration(spec, copy(Xoshiro(s)))

				#	Run the simulator with the same generator state
					res = sirdif(net.alst, net.vlst, [1], 0.02, spec, maxtime, 0.5,
								 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
								 rng = Xoshiro(s))
					ilog = res["infection_log"]

				#	Compare the observed recovery day with the prediction
					if D <= maxtime
						infected_rows = ilog[ilog.time .< D, :n_infected]
						recovery_row = ilog[ilog.time .== D, :]
						ok = all(infected_rows .== 1) && nrow(recovery_row) == 1 &&
							 recovery_row.n_infected[1] == 0 && recovery_row.n_recovered[1] == 1
					else
						ok = all(ilog.n_infected .== 1)
					end

				#	Record the case
					results.scale[k] = scale
					results.shape[k] = shape
					results.seed[k] = s
					results.drawn_days[k] = D
					results.passed[k] = ok
			end

		#	Report
			all_passed = all(results.passed)
			if verbose
				println("=" ^ 60)
				println("Weibull Recovery Timing (isolated node)")
				println("=" ^ 60)
				println("  Cases: $(n_cases)   Passed: $(sum(results.passed))")
				println("  Shortest draw: $(minimum(results.drawn_days))   Longest draw: $(maximum(results.drawn_days))")
				println("RESULT: $(all_passed ? "ALL CASES PASSED" : "SOME CASES FAILED")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed)
	end
	@doc raw"""
	**Description**
	Verify end to end that a Weibull recovery draw sets the day an infected individual recovers. A single isolated seed makes the recovery draw the first use of the generator, so the expected duration can be computed from a copy of it.

	**Usage**
	`test_weibull_recovery_timing(; scales=[0.5, 10.0], shapes=[0.8, 2.0], seeds=1:50, maxtime=400, verbose=true)`

	**Value**
	A `NamedTuple` with the per-case `results::DataFrame` and `all_passed::Bool`.

	**See Also**
	`sirdif`, `draw_duration`, `WeibullDuration`
	""" test_weibull_recovery_timing

#	Recovery and Immunity Schedule on a Deterministic Path Network
	function test_waning_schedule(; L::Int = 10, D::Int = 3, M::Int = 2, maxtime::Int = 50,
									verbose::Bool = true)
		"""
		Args:
			L::Int: number of nodes on the path (default = 10)
			D::Int: fixed recovery duration in days (default = 3)
			M::Int: fixed immunity duration in days (default = 2)
			maxtime::Int: simulation horizon (default = 50)
			verbose::Bool: print a report (default = true)
		Returns:
			NamedTuple: (checks::Dict{String,Bool}, all_passed::Bool)
		Notes:
			Node 1 is an isolated seed. Nodes 2 to L+1 form a path seeded at
			node 2. With b_int = 50 and inf_r = 1, every contact is made and
			every contact transmits, so node k+1 on the path is infected on
			day k-1. The expected counts follow exactly from D and M, which
			checks the recovery countdown, the immunity schedule (M full
			immune days), and the cumulative incidence column. Both call
			forms are run and must agree.
		"""

		#	Build the network: isolated node 1, path 2 - 3 - ... - (L+1)
			neighbors = Vector{Vector{Int}}(undef, L + 1)
			neighbors[1] = Int[]
			for k in 1:L
				id = k + 1
				nbrs = Int[]
				if k > 1
					push!(nbrs, id - 1)
				end
				if k < L
					push!(nbrs, id + 1)
				end
				neighbors[id] = nbrs
			end
			net = build_test_network(neighbors)
			n = L + 1

		#	Expected schedule
			infection_days = vcat([0], collect(0:(L - 1)))
			t_end = (L - 1) + D
			expected = expected_path_counts(infection_days, D, M, t_end)

		#	Run both call forms with certain contact and transmission
			res_spec = sirdif(net.alst, net.vlst, [1, 2], 1.0, FixedDuration(D), maxtime, 0.5,
							  50.0, 0.0, 0.0, 0.0, 0.0, 0.0;
							  immunity = FixedDuration(M), rng = Xoshiro(7))
			res_int = sirdif(net.alst, net.vlst, [1, 2], 1.0, D, maxtime, 0.5,
							 50.0, 0.0, 0.0, 0.0, 0.0, 0.0;
							 immunity_duration = M, rng = Xoshiro(7))
			ilog = res_spec["infection_log"]

		#	Compare observed and expected counts for days 1 to t_end
			checks = Dict{String, Bool}()
			obs = ilog[(ilog.time .>= 1) .& (ilog.time .<= t_end), :]
			exp_rows = expected[expected.time .>= 1, :]
			checks["row_count"] = nrow(obs) == nrow(exp_rows)
			if checks["row_count"]
				checks["n_infected"] = all(obs.n_infected .== exp_rows.n_infected)
				checks["n_recovered"] = all(obs.n_recovered .== exp_rows.n_recovered)
				checks["cumulative"] = all(isapprox.(obs.prop_cum_infected, exp_rows.n_ever ./ n; atol = 1e-12))
			end

		#	Check the initial row, the extinction row, and the call forms
			checks["initial_row"] = ilog.n_infected[1] == 2 && isapprox(ilog.prop_cum_infected[1], 2 / n; atol = 1e-12)
			checks["extinction_row"] = ilog.time[end] == maxtime && ilog.n_infected[end] == 0
			checks["prop_ever_declines"] = any(diff(ilog.prop_ever_infected[ilog.time .<= t_end]) .< 0)
			checks["call_forms_agree"] = isequal(res_spec["infection_log"], res_int["infection_log"]) &&
										 res_spec["final_state"] == res_int["final_state"]

		#	Report
			all_passed = all(values(checks))
			if verbose
				println("=" ^ 60)
				println("Recovery and Immunity Schedule (path network, D=$(D), M=$(M))")
				println("=" ^ 60)
				for (name, passed) in sort(collect(checks))
					println("  $(rpad(name, 22)) $(passed ? "PASS" : "FAIL")")
				end
				println("RESULT: $(all_passed ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")")
			end

		#	Assembling Result
			return (checks = checks, all_passed = all_passed, ilog = ilog, expected = expected)
	end
	@doc raw"""
	**Description**
	Verify the recovery countdown, the waning immunity schedule, and the cumulative incidence column against exact expected counts on a deterministic path network.

	**Usage**
	`test_waning_schedule(; L=10, D=3, M=2, maxtime=50, verbose=true)`

	**Details**
	With `b_int = 50` and `inf_r = 1`, the contact probability rounds to 1.0 and the transmission probability is 1, so the epidemic moves one step along the path each day regardless of the random draws. An individual who recovers at the end of day r must be counted as recovered in rows r to r+M-1 and as susceptible from row r+M. The validated simulator's fixed-immunity path gave M-1 immune days; this test confirms the correction.

	**Value**
	A `NamedTuple` with `checks`, `all_passed`, the observed `log`, and the `expected` counts.

	**See Also**
	`sirdif`, `expected_path_counts`
	""" test_waning_schedule

#	UT-1.6: Simulator Invariants
	function test_simulator_invariants(network_data::NamedTuple; n_runs::Int = 100, maxtime::Int = 200,
									   master_seed::Int = 2026, verbose::Bool = true)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
			n_runs::Int: number of runs (default = 100)
			maxtime::Int: simulation horizon (default = 200)
			master_seed::Int: run r uses Xoshiro(master_seed + r) (default = 2026)
			verbose::Bool: print a report (default = true)
		Returns:
			NamedTuple: (results::DataFrame, all_passed::Bool)
		Notes:
			Weibull recovery (scale 10, shape 2; every tenth run uses scale 0.5
			to exercise the floor) and Weibull immunity (scale 30, shape 2).
			One random seed node per run. Behavioral parameters follow the R
			settings in test_sir_model. Criterion (e), daily expansion by
			summarize.jl, is added when summarize.jl exists.
		"""

		#	Network
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)
			ids = alst[:, 1]

		#	Pre-allocate results
			results = DataFrame(run = collect(1:n_runs), rec_scale = zeros(n_runs),
								cum_monotone = falses(n_runs), cum_bounded = falses(n_runs),
								cum_dominates = falses(n_runs), counts_valid = falses(n_runs),
								counters_valid = falses(n_runs), states_valid = falses(n_runs),
								prop_ever_declined = falses(n_runs), final_cum = zeros(n_runs))

		#	Iterate over runs
			for r in 1:n_runs
				#	Draw the seed node and set durations
					rng = Xoshiro(master_seed + r)
					seed_id = ids[rand(rng, 1:n)]
					rec_scale = (r % 10 == 0) ? 0.5 : 10.0
					recovery = WeibullDuration(rec_scale, 2.0)
					immunity = WeibullDuration(30.0, 2.0)

				#	Run the simulator
					res = sirdif(alst, vlst, [seed_id], 0.33, recovery, maxtime, 0.5,
								 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
								 immunity = immunity, rng = rng)
					ilog = res["infection_log"]
					st = res["final_state"]

				#	(a) Cumulative incidence is non-decreasing, at most 1, and at least (NI + NR)/N
					cum = ilog.prop_cum_infected
					results.cum_monotone[r] = all(diff(cum) .>= 0)
					results.cum_bounded[r] = all(cum .<= 1.0 + 1e-12)
					results.cum_dominates[r] = all(cum .>= ilog.prop_ever_infected .- 1e-12)

				#	Counts are within bounds and time increases
					results.counts_valid[r] = all(ilog.n_infected .>= 0) && all(ilog.n_recovered .>= 0) &&
											  all(ilog.n_infected .+ ilog.n_recovered .<= n) &&
											  all(diff(ilog.time) .> 0)

				#	(b, c) Counters: infected have 1 or more days left, recovered have 1 or more immune days left, susceptible have none
					infected = st[:, 2] .== 1
					recovered = st[:, 4] .== 1
					susceptible = st[:, 3] .== 1
					results.counters_valid[r] = all(st[infected, 5] .>= 1) && all(st[recovered, 8] .>= 1) &&
												all(st[susceptible, 5] .== 0) && all(st[:, 5] .>= 0) &&
												all(st[:, 8] .>= 0)

				#	(d) Each individual occupies exactly one state
					results.states_valid[r] = all(st[:, 2] .+ st[:, 3] .+ st[:, 4] .== 1)

				#	Record whether waning was active and the final size
					results.prop_ever_declined[r] = any(diff(ilog.prop_ever_infected) .< 0)
					results.rec_scale[r] = rec_scale
					results.final_cum[r] = cum[end]
			end

		#	Report
			criteria = [:cum_monotone, :cum_bounded, :cum_dominates, :counts_valid, :counters_valid, :states_valid]
			all_passed = all(all(results[!, c]) for c in criteria)
			if verbose
				println("=" ^ 60)
				println("UT-1.6: Simulator Invariants ($(n_runs) runs, $(n) nodes)")
				println("=" ^ 60)
				for c in criteria
					println("  $(rpad(String(c), 18)) $(sum(results[!, c])) / $(n_runs)")
				end
				println("  Runs where (NI + NR)/N declined (waning active): $(sum(results.prop_ever_declined))")
				println("  Median final cumulative incidence: $(round(median(results.final_cum), digits=4))")
				println("  Criterion (e) pending summarize.jl")
				println("RESULT: $(all_passed ? "ALL RUNS PASSED" : "SOME RUNS FAILED")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed)
	end
	@doc raw"""
	**Description**
	UT-1.6. Run the simulator repeatedly with Weibull recovery and waning immunity and check the invariants listed in the roadmap (§4.8).

	**Usage**
	`test_simulator_invariants(network_data; n_runs=100, maxtime=200, master_seed=2026, verbose=true)`

	**Details**
	Per run: (a) `prop_cum_infected` is non-decreasing, at most 1, and never below `prop_ever_infected`; counts stay within bounds and time increases; (b, c) no negative counters, every infected individual has at least one day to recovery, and every recovered individual has at least one immune day left; (d) each individual occupies exactly one of S, I, and R. Criterion (e) is added with `summarize.jl`. The report also counts runs in which (NI + NR)/N declined, which confirms that waning immunity was active.

	**Value**
	A `NamedTuple` with the per-run `results::DataFrame` and `all_passed::Bool`.

	**See Also**
	`sirdif`, `test_waning_schedule`, `test_weibull_recovery_timing`
	""" test_simulator_invariants

#	Helper Function for test_regression: compare a new result with the reference result
	function compare_to_reference(new_res::Dict, ref_res::Dict)
		"""
		Args:
			new_res::Dict: result from diffusion_sim.sirdif
			ref_res::Dict: result from ReferenceSirdif.sirdif
		Returns:
			NamedTuple: (log_match::Bool, state_match::Bool, cum_match::Bool)
		Notes:
			The new log must reproduce the reference's six columns exactly and
			append prop_cum_infected equal to prop_ever_infected. The new final
			state must reproduce the reference's seven columns exactly.
		"""

		#	Compare logs column by column
			new_log = new_res["infection_log"]
			ref_log = ref_res["infection_log"]
			ref_cols = names(ref_log)
			log_match = names(new_log)[1:length(ref_cols)] == ref_cols &&
						nrow(new_log) == nrow(ref_log) &&
						all(isequal(new_log[!, c], ref_log[!, c]) for c in ref_cols)

		#	Compare final states on the reference columns
			new_state = new_res["final_state"]
			ref_state = ref_res["final_state"]
			state_match = new_state[:, 1:size(ref_state, 2)] == ref_state

		#	Appended column equals prop_ever_infected under permanent immunity
			cum_match = names(new_log)[end] == "prop_cum_infected" &&
						isequal(new_log.prop_cum_infected, new_log.prop_ever_infected)

		#	Assembling Result
			return (log_match = log_match, state_match = state_match, cum_match = cum_match)
	end

#	UT-1.7: Regression Against the Validated Simulator
	function test_regression(network_data::NamedTuple, network_sas::NamedTuple;
							 n_iter::Int = 10, verbose::Bool = true)
		"""
		Args:
			network_data::NamedTuple: sim_prep output without SAS preprocessing (test_sir_model setting)
			network_sas::NamedTuple: sim_prep output with SAS preprocessing (replication setting)
			n_iter::Int: iterations per coefficient setting in the replication part (default = 10)
			verbose::Bool: print a report (default = true)
		Returns:
			NamedTuple: (results::DataFrame, all_passed::Bool, stream_equivalence::Bool)
		Notes:
			Constant recovery and permanent immunity throughout. Each pair of
			runs reseeds the global generator identically before calling the
			reference and the new simulator. Part C checks that the
			DurationSpec call form matches the integer form. The last check,
			whether Xoshiro(s) reproduces Random.seed!(s), is informative and
			not part of the pass criterion.
		"""

		#	Pre-allocate results
			results = DataFrame(part = String[], case = String[], log_match = Bool[],
								state_match = Bool[], cum_match = Bool[])

		#	Part A: test_sir_model settings (five fixed seed nodes)
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)
			seed_nodes = [s for s in [5, 10, 15, 20, 25] if 1 <= s <= n]
			for s in [42, 123, 2026]
				#	Reference run
					Random.seed!(s)
					ref_res = ReferenceSirdif.sirdif(alst, vlst, seed_nodes, 0.33, 14, 200, 0.5,
												 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
												 transmission_method = :weighted, immunity_duration = nothing)
				#	New run
					Random.seed!(s)
					new_res = sirdif(alst, vlst, seed_nodes, 0.33, 14, 200, 0.5,
								 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
								 transmission_method = :weighted, immunity_duration = nothing)
				#	Record
					comparison = compare_to_reference(new_res, ref_res)
					push!(results, ("A: test_sir_model", "seed $(s)", comparison.log_match, comparison.state_match, comparison.cum_match))
			end

		#	Part B: replication settings (one random seed node per run, as in replicate_sas_simulation)
			alst_s = network_sas.alst
			vlst_s = network_sas.vlst
			n_s = size(alst_s, 1)
			design = [(0.0, 0.0, 0.0), (-0.5, -3.5, -1.5), (-1.0, -6.5, -3.5)]
			for (i, (b_peer, b_glob, b_self)) in enumerate(design)
				for it in 1:n_iter
					#	Reference run
						Random.seed!(10_000 * i + it)
						seednode = alst_s[rand(1:n_s), 1]
						ref_res = ReferenceSirdif.sirdif(alst_s, vlst_s, [seednode], 0.02, 14, 200, 0.75,
													 -0.1, 1.0, b_peer, b_glob, b_self, -0.1;
													 transmission_method = :weighted, immunity_duration = nothing)
					#	New run
						Random.seed!(10_000 * i + it)
						seednode_new = alst_s[rand(1:n_s), 1]
						new_res = sirdif(alst_s, vlst_s, [seednode_new], 0.02, 14, 200, 0.75,
									 -0.1, 1.0, b_peer, b_glob, b_self, -0.1;
									 transmission_method = :weighted, immunity_duration = nothing)
					#	Record
						comparison = compare_to_reference(new_res, ref_res)
						push!(results, ("B: replication", "design $(i), iter $(it)",
										comparison.log_match && seednode == seednode_new, comparison.state_match, comparison.cum_match))
				end
			end

		#	Part C: DurationSpec call form matches the integer form
			Random.seed!(42)
			int_form = sirdif(alst, vlst, seed_nodes, 0.33, 14, 200, 0.5,
							  -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
			Random.seed!(42)
			spec_form = sirdif(alst, vlst, seed_nodes, 0.33, FixedDuration(14), 200, 0.5,
							   -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
			forms_match = isequal(int_form["infection_log"], spec_form["infection_log"]) &&
						  int_form["final_state"] == spec_form["final_state"]
			push!(results, ("C: call forms", "seed 42", forms_match, forms_match, forms_match))

		#	Informative: explicit Xoshiro(s) versus Random.seed!(s)
			Random.seed!(42)
			global_form = sirdif(alst, vlst, seed_nodes, 0.33, 14, 200, 0.5,
								 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
			explicit_form = sirdif(alst, vlst, seed_nodes, 0.33, 14, 200, 0.5,
								   -0.1, 1.0, -0.5, -3.5, -1.5, -0.1; rng = Xoshiro(42))
			stream_equivalence = isequal(global_form["infection_log"], explicit_form["infection_log"])

		#	Report
			all_passed = all(results.log_match) && all(results.state_match) && all(results.cum_match)
			if verbose
				println("=" ^ 60)
				println("UT-1.7: Regression Against the Validated Simulator")
				println("=" ^ 60)
				for g in groupby(results, :part)
					n_ok = sum(g.log_match .& g.state_match .& g.cum_match)
					println("  $(rpad(first(g.part), 20)) $(n_ok) / $(nrow(g)) cases identical")
				end
				println("  Informative: Xoshiro(42) reproduces Random.seed!(42): $(stream_equivalence ? "YES" : "NO")")
				println("RESULT: $(all_passed ? "ALL CASES IDENTICAL" : "SOME CASES DIFFER")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed, stream_equivalence = stream_equivalence)
	end
	@doc raw"""
	**Description**
	UT-1.7. Run the frozen reference copy of the validated `sirdif` and the modified `sirdif` side by side with identical seeds, under constant recovery and permanent immunity, and require identical output.

	**Usage**
	`test_regression(network_data, network_sas; n_iter=10, verbose=true)`

	**Details**
	Part A repeats the `test_sir_model` settings with three global seeds. Part B repeats the `replicate_sas_simulation` pattern (one random seed node per run, drawn from the global generator) on the SAS-preprocessed network for three coefficient settings. Part C checks that `FixedDuration(14)` and the integer form `14` give identical results. For every case, the six reference log columns and the seven reference state columns must match exactly, and the appended `prop_cum_infected` must equal `prop_ever_infected`.

	**Value**
	A `NamedTuple` with the per-case `results::DataFrame`, `all_passed::Bool`, and the informative `stream_equivalence::Bool`.

	**See Also**
	`ReferenceSirdif.sirdif`, `sirdif`, `compare_to_reference`
	""" test_regression

#############
#   TESTS   #
#############

#	Network Paths (EC2 test data)
	graphml_path = "/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml"

#	Load Networks
	network_data = sim_prep(graphml_path)
	network_sas = sim_prep(graphml_path; sas_transformation = true)

#	Duration Draws and Deterministic Schedules
	draw_report = test_draw_duration()
	timing_report = test_weibull_recovery_timing()
	waning_report = test_waning_schedule()

#	UT-1.6: Simulator Invariants
	invariant_report = test_simulator_invariants(network_data)

#	UT-1.7: Regression Against the Validated Simulator
	regression_report = test_regression(network_data, network_sas)

#	Summary
	println("=" ^ 60)
	println("Simulator Test Summary")
	println("=" ^ 60)
	println("  Duration draws:         $(draw_report.all_passed ? "PASS" : "FAIL")")
	println("  Weibull recovery timing: $(timing_report.all_passed ? "PASS" : "FAIL")")
	println("  Waning schedule:        $(waning_report.all_passed ? "PASS" : "FAIL")")
	println("  UT-1.6 invariants:      $(invariant_report.all_passed ? "PASS" : "FAIL")")
	println("  UT-1.7 regression:      $(regression_report.all_passed ? "PASS" : "FAIL")")

#Run Summary Tests: Daily Expansion, Summaries, Classification, and UT-1.6(e)
#Jonathan H. Morgan, Ph.D.
#12 October 2026
#Revised 12 October 2026: t_i is the inflection time (the day of peak incidence),
#the t_peak field is removed, and fit_rmse is checked. The shape a now runs from
#-1 to 10; frac_at_inflection is checked, and the curve diagnostics add the start
#time, a smoothed incidence peak, and indicators of persistent infection.

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

#   Checking for LsqFit (required by richards.jl)
	if Base.find_package("LsqFit") === nothing
		error("LsqFit is not in this environment. Run: using Pkg; Pkg.add(\"LsqFit\")")
	end

################
#   PACKAGES   #
################

#	Loading Packages
	using DataFrames
	using Random
	using Statistics
	using diffusion_sim

#   Loading Test Modules: Richards and Summary
	include("/workspace/ergm-based-diffusion-simulations/julia_env/src/richards.jl")
	include("/workspace/ergm-based-diffusion-simulations/julia_env/src/summarize.jl")
	using .richards_model
	using .run_summary

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Summary Tests: report a dictionary of checks
	function report_checks(title::String, checks::Dict{String, Bool})
		"""
		Args:
			title::String: report heading
			checks::Dict{String, Bool}: named pass/fail results
		Returns:
			Bool: true if every check passed
		"""

		#	Print each check in name order
			println("=" ^ 60)
			println(title)
			println("=" ^ 60)
			for (name, passed) in sort(collect(checks))
				println("  $(rpad(name, 40)) $(passed ? "PASS" : "FAIL")")
			end
			all_passed = all(values(checks))
			println("RESULT: $(all_passed ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")")

		#	Assembling Result
			return all_passed
	end

#	Helper Function for the Summary Tests: check that a call throws a given error type
	function throws_error(f::Function, err_type::Type)
		"""
		Args:
			f::Function: zero-argument function to call
			err_type::Type: expected exception type
		Returns:
			Bool: true if f() throws an exception of err_type
		"""

		#	Call and inspect
			try
				f()
				return false
			catch e
				return e isa err_type
			end
	end

#	Helper Function for the Summary Tests: build a small network in sim_prep format
	function build_test_network(neighbors::Vector{Vector{Int}})
		"""
		Args:
			neighbors::Vector{Vector{Int}}: neighbors[i] lists the IDs tied to node i (IDs are 1:n)
		Returns:
			NamedTuple: (alst::Matrix{Int}, vlst::Matrix{Float64})
		Notes:
			Same construction as in simulator_tests.jl. Ties must be listed in
			both directions; all weights are 1.0.
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

#	Helper Function for the Summary Tests: synthetic infection log
	function synthetic_log(times::Vector{Float64}, n_infected::Vector{Float64}, cum::Vector{Float64})
		"""
		Args:
			times::Vector{Float64}: day of each row
			n_infected::Vector{Float64}: infected count on each row
			cum::Vector{Float64}: cumulative proportion on each row
		Returns:
			DataFrame: a log with the seven sirdif columns
		Notes:
			Columns not used by the summaries are filled with zeros.
		"""

		#	Assemble the columns in sirdif order
			n = length(times)
			return DataFrame(time = times, n_infected = n_infected, prop_ever_infected = cum,
							 prop_currently_infected = zeros(n), n_recovered = zeros(n),
							 prop_recovered = zeros(n), prop_cum_infected = cum)
	end

#	expand_daily_log: Filling the Days After Extinction
	function test_expand_daily_log()
		"""
		Returns:
			Bool: true if every check passed
		Notes:
			Uses synthetic logs: a complete log, an extinction log with the
			single SAS-style row at maxtime, and four invalid logs.
		"""

		#	Initialize checks
			checks = Dict{String, Bool}()

		#	Complete log is unchanged
			full = synthetic_log(collect(0.0:5.0), [1.0, 2, 3, 2, 1, 1], [0.1, 0.2, 0.3, 0.3, 0.3, 0.3])
			daily_full = expand_daily_log(full, 5)
			checks["complete log unchanged"] = isequal(daily_full, full)

		#	Extinction log: days 0 to 3, then one row at maxtime = 10
			ext = synthetic_log([0.0, 1.0, 2.0, 3.0, 10.0], [1.0, 2.0, 1.0, 0.0, 0.0],
								[0.1, 0.2, 0.3, 0.3, 0.3])
			daily_ext = expand_daily_log(ext, 10)
			checks["extinction log has maxtime + 1 rows"] = nrow(daily_ext) == 11
			checks["days are 0 to maxtime"] = daily_ext.time == collect(0.0:10.0)
			checks["filled days carry the last state"] = all(daily_ext.n_infected[5:10] .== 0.0) &&
														 all(daily_ext.prop_cum_infected[5:11] .== 0.3)
			checks["simulated days unchanged"] = daily_ext.n_infected[1:4] == [1.0, 2.0, 1.0, 0.0]
			checks["input log not modified"] = nrow(ext) == 5

		#	Invalid logs
			checks["empty log rejected"] = throws_error(() -> expand_daily_log(synthetic_log(Float64[], Float64[], Float64[]), 5), ArgumentError)
			checks["log not starting at 0 rejected"] = throws_error(() -> expand_daily_log(synthetic_log([1.0, 2.0], [1.0, 1.0], [0.1, 0.1]), 5), ArgumentError)
			checks["fractional day rejected"] = throws_error(() -> expand_daily_log(synthetic_log([0.0, 1.5], [1.0, 1.0], [0.1, 0.1]), 5), ArgumentError)
			checks["log past maxtime rejected"] = throws_error(() -> expand_daily_log(synthetic_log([0.0, 9.0], [1.0, 1.0], [0.1, 0.1]), 5), ArgumentError)

		#	Assembling Result
			return report_checks("expand_daily_log: Filling the Days After Extinction", checks)
	end

#	trajectory_summaries: Hand-Computed Values
	function test_trajectory_summaries()
		"""
		Returns:
			Bool: true if every check passed
		"""

		#	A short series with a known peak, slope, and end
			I = [0.01, 0.03, 0.08, 0.10, 0.06, 0.02, 0.0, 0.0]
			C = [0.01, 0.03, 0.09, 0.16, 0.18, 0.19, 0.19, 0.19]
			s = trajectory_summaries(I, C; surge_time=2)
			checks = Dict{String, Bool}()
			checks["peak infected"] = s.peak_infected == 0.10
			checks["peak day"] = s.peak_day == 3
			checks["infected at surge day"] = s.infected_at_surge == 0.08
			checks["max slope"] = isapprox(s.max_slope, 0.05; atol=1e-15)
			checks["final size"] = s.final_size == 0.19
			checks["duration"] = s.duration == 6 && s.extinguished

		#	Persistent infection and a surge beyond the series
			s2 = trajectory_summaries([0.01, 0.02, 0.03], [0.01, 0.02, 0.04]; surge_time=21)
			checks["censored duration"] = s2.duration == 2 && !s2.extinguished
			checks["surge beyond series is NaN"] = isnan(s2.infected_at_surge)
			checks["unequal lengths rejected"] = throws_error(() -> trajectory_summaries([0.1], [0.1, 0.2]), ArgumentError)

		#	Assembling Result
			return report_checks("trajectory_summaries: Hand-Computed Values", checks)
	end

#	summarize_run: Deterministic Networks for Each Run Class
	function test_run_classes()
		"""
		Returns:
			Bool: true if every check passed
		Notes:
			With b_int = 50 and inf_r = 1, every contact is made and transmits,
			so these runs are deterministic.
			Non-starter: two isolated nodes, seed node 1.
			Extinct: 100 nodes; the seed's only tie is node 2, which has no
			other ties, so 2 of 100 are infected (growth 0.01 beyond C0).
			Path: an isolated seed plus a 10-node path seeded at one end,
			recovery 3 days, immunity 2 days (the waning-schedule network
			from simulator_tests.jl), with exact expected summaries.
		"""

		#	Common settings
			checks = Dict{String, Bool}()
			region = TargetRegion(K_max=1.0)
			maxtime = 50

		#	Non-starter
			net_ns = build_test_network([Int[], Int[]])
			res_ns = sirdif(net_ns.alst, net_ns.vlst, [1], 1.0, FixedDuration(3), maxtime, 0.5,
							50.0, 0.0, 0.0, 0.0, 0.0, 0.0; rng=Xoshiro(1))
			for policy in (:exclude, :success, :failure)
				s = summarize_run(res_ns["infection_log"], 2, maxtime, region; extinct_policy=policy)
				checks["non-starter ($(policy))"] = s.class == :nonstarter && !s.viable &&
													s.in_target === nothing && s.C0 == 0.5 && isnan(s.fit_rmse)
			end

		#	Extinct
			nbrs = [Int[] for _ in 1:100]
			nbrs[1] = [2]
			nbrs[2] = [1]
			net_ext = build_test_network(nbrs)
			res_ext = sirdif(net_ext.alst, net_ext.vlst, [1], 1.0, FixedDuration(3), maxtime, 0.5,
							 50.0, 0.0, 0.0, 0.0, 0.0, 0.0; rng=Xoshiro(2))
			ilog_ext = res_ext["infection_log"]
			s_ex = summarize_run(ilog_ext, 100, maxtime, region; extinct_policy=:exclude)
			s_fa = summarize_run(ilog_ext, 100, maxtime, region; extinct_policy=:failure)
			s_su = summarize_run(ilog_ext, 100, maxtime, region; extinct_policy=:success)
			checks["extinct: class and final size"] = s_ex.class == :extinct && s_ex.viable &&
													  isapprox(s_ex.final_size, 0.02; atol=1e-15) && isnan(s_ex.fit_rmse)
			checks["extinct: :exclude gives nothing"] = s_ex.in_target === nothing
			checks["extinct: :failure gives false"] = s_fa.in_target === false
			checks["extinct: :success within K_max"] = s_su.in_target === true

		#	Path network with waning immunity
			L = 10
			nbrs_p = Vector{Vector{Int}}(undef, L + 1)
			nbrs_p[1] = Int[]
			for k in 1:L
				id = k + 1
				nbrs_p[id] = vcat(k > 1 ? [id - 1] : Int[], k < L ? [id + 1] : Int[])
			end
			net_p = build_test_network(nbrs_p)
			res_p = sirdif(net_p.alst, net_p.vlst, [1, 2], 1.0, FixedDuration(3), maxtime, 0.5,
						   50.0, 0.0, 0.0, 0.0, 0.0, 0.0; immunity=FixedDuration(2), rng=Xoshiro(3))
			s_p = summarize_run(res_p["infection_log"], 11, maxtime, region)
			checks["path: C0 = 2/11"] = s_p.C0 == 2 / 11
			checks["path: peak 4/11 on day 2"] = s_p.peak_infected == 4 / 11 && s_p.peak_day == 2
			checks["path: max slope 1/11"] = isapprox(s_p.max_slope, 1 / 11; atol=1e-15)
			checks["path: final size 1"] = s_p.final_size == 1.0
			checks["path: ends on day 12"] = s_p.duration == 12 && s_p.extinguished
			checks["path: maxtime + 1 days"] = s_p.n_days == maxtime + 1
			checks["path: class matches status"] = (s_p.status == :ok) == (s_p.class == :established)
			checks["record has fit_rmse, no t_peak"] = haskey(s_p, :fit_rmse) && !haskey(s_p, :t_peak)
			checks["record has frac_at_inflection"] = haskey(s_p, :frac_at_inflection)

		#	Report the path fit status for the record
			println("  (path network fit status: $(s_p.status))")

		#	Assembling Result
			return report_checks("summarize_run: Deterministic Run Classes", checks)
	end

#	UT-1.6(e) and Pipeline Consistency on a Village Network
	function test_pipeline_on_network(network_data::NamedTuple; n_runs::Int=100, maxtime::Int=200,
									  master_seed::Int=2026, verbose::Bool=true)
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
			Same settings as UT-1.6 (Weibull recovery scale 10, shape 2, with
			scale 0.5 on every tenth run; Weibull immunity scale 30, shape 2;
			R-setting coefficients), with one seed drawn from nodes that have
			at least one tie. The region is permissive (K_max = 1), so every
			established run must be in the target.
			Per run: the expanded log has maxtime + 1 consecutive days
			(UT-1.6 criterion e), C0 = 1/N, the final size matches the log,
			and the class, status, viability, and in_target agree. For an
			established run the inflection time t_i lies within [0, maxtime],
			a lies within [-1, 10], frac_at_inflection lies within [0, 1),
			and fit_rmse is finite; a non-starter or extinct run has
			fit_rmse = NaN.
		"""

		#	Network and eligible seeds
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)
			ids = alst[:, 1]
			eligible = findall(i -> alst[i, 2] > 0, 1:n)
			region = TargetRegion(K_max=1.0)

		#	Pre-allocate results
			results = DataFrame(run = Int[], rec_scale = Float64[], status = Symbol[], class = Symbol[],
								final_size = Float64[], K = Float64[], t_i = Float64[], a = Float64[],
								frac = Float64[], fit_rmse = Float64[], in_target = Any[],
								days_ok = Bool[], c0_ok = Bool[], final_ok = Bool[], class_ok = Bool[],
								curve_ok = Bool[])

		#	Iterate over runs
			for r in 1:n_runs
				#	Simulate
					rng = Xoshiro(master_seed + r)
					seed_id = ids[eligible[rand(rng, 1:length(eligible))]]
					rec_scale = (r % 10 == 0) ? 0.5 : 10.0
					res = sirdif(alst, vlst, [seed_id], 0.33, WeibullDuration(rec_scale, 2.0), maxtime, 0.5,
								 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
								 immunity=WeibullDuration(30.0, 2.0), rng=rng)
					ilog = res["infection_log"]

				#	Summarize
					s = summarize_run(ilog, n, maxtime, region; rec_scale=rec_scale, rec_shape=2.0)
					daily = expand_daily_log(ilog, maxtime)

				#	(e) one row per day, consecutive
					days_ok = nrow(daily) == maxtime + 1 && daily.time == collect(0.0:maxtime)

				#	C0 and final size
					c0_ok = s.C0 == 1 / n
					final_ok = s.final_size == ilog.prop_cum_infected[end] && s.final_size >= s.C0

				#	Class, status, viability, and in_target agree
					class_ok = if s.class == :nonstarter
						s.status == :nonstarter && !s.viable && s.in_target === nothing && s.final_size == s.C0
					elseif s.class == :extinct
						s.status == :extinct && s.viable && s.in_target === nothing
					elseif s.class == :established
						s.status == :ok && s.viable && s.in_target === true
					else
						!s.fit_ok && s.in_target === nothing
					end

				#	Inflection within the horizon and residual reported
					curve_ok = if s.class == :established
						0.0 <= s.t_i <= maxtime && -1.0 <= s.a <= 10.0 &&
						0.0 <= s.frac_at_inflection < 1.0 && isfinite(s.fit_rmse)
					elseif s.class == :nonstarter || s.class == :extinct
						isnan(s.fit_rmse)
					else
						true
					end

				#	Record
					push!(results, (r, rec_scale, s.status, s.class, s.final_size, s.K, s.t_i, s.a,
									s.frac_at_inflection, s.fit_rmse, s.in_target, days_ok, c0_ok, final_ok, class_ok, curve_ok))
			end

		#	Report
			criteria = [:days_ok, :c0_ok, :final_ok, :class_ok, :curve_ok]
			all_passed = all(all(results[!, c]) for c in criteria)
			if verbose
				println("=" ^ 60)
				println("UT-1.6(e) and Pipeline Consistency ($(n_runs) runs, $(n) nodes)")
				println("=" ^ 60)
				for c in criteria
					println("  $(rpad(String(c), 12)) $(sum(results[!, c])) / $(n_runs)")
				end
				println("  Class counts:")
				for g in groupby(results, :class)
					println("    $(rpad(String(first(g.class)), 12)) $(nrow(g))")
				end
				println("  Status counts:")
				for g in groupby(results, :status)
					println("    $(rpad(String(first(g.status)), 16)) $(nrow(g))")
				end
				n_viable = count(results.class .!= :nonstarter)
				n_est = count(results.class .== :established)
				n_fail = count(results.class .== :fit_failed)
				println("  Established share of viable runs: $(n_est) / $(n_viable)")
				println("  Fit success among runs that were fitted: $(n_est) / $(n_est + n_fail)")
				est = results[results.class .== :established, :]
				if nrow(est) > 0
					println("  Established runs: median K = $(round(median(est.K), digits=4)), " *
							"median final size = $(round(median(est.final_size), digits=4))")
					println("  Established runs: median t_i = $(round(median(est.t_i), digits=2)), " *
							"median fit_rmse = $(round(median(est.fit_rmse), sigdigits=3)), " *
							"a < 0 in $(count(est.a .< 0.0)) / $(nrow(est)), a = -1 in $(count(est.a .<= -1.0 + 1.0e-6))")
					println("  Established runs: median a = $(round(median(est.a), digits=3)), " *
							"median frac_at_inflection = $(round(median(est.frac), digits=3))")
				end
				println("RESULT: $(all_passed ? "ALL RUNS CONSISTENT" : "SOME RUNS INCONSISTENT")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed)
	end

#############
#   TESTS   #
#############

#	Network Path (EC2 test data)
	graphml_path = "/workspace/ergm_data/sim_test_data/WeakCore1_3.2_2.graphml"
	network_data = sim_prep(graphml_path)

#	Run Each Test
	expand_passed = test_expand_daily_log()
	traj_passed = test_trajectory_summaries()
	classes_passed = test_run_classes()
	pipeline_report = test_pipeline_on_network(network_data)

#	Summary
	println("=" ^ 60)
	println("Run Summary Test Summary")
	println("=" ^ 60)
	println("  expand_daily_log:         $(expand_passed ? "PASS" : "FAIL")")
	println("  trajectory_summaries:     $(traj_passed ? "PASS" : "FAIL")")
	println("  run classes:              $(classes_passed ? "PASS" : "FAIL")")
	println("  UT-1.6(e) and pipeline:   $(pipeline_report.all_passed ? "PASS" : "FAIL")")

###################
#   CURVE TESTS   #
###################

#	Testing Eligibility Rules: Fitted Inflection Against the Observed Incidence Peak
#	obs_inc_peak is the midpoint of the largest one-day rise in cumulative incidence;
#	obs_inc_peak_ma7 is the same peak after a 7-day centered moving average;
#	obs_peak_day is the day of peak prevalence, which follows the incidence peak;
#	t_start is the fitted start, richards_start (-Inf for a ≥ 0);
#	extinguished and late_prev (mean prevalence on days 150 to 200) flag persistent infection
	alst = network_data.alst; vlst = network_data.vlst; n = size(alst, 1)
	eligible = findall(i -> alst[i, 2] > 0, 1:n); ids = alst[:, 1]
	rows = []
	for r in 1:30
		rng = Xoshiro(2026 + r)
		seed_id = ids[eligible[rand(rng, 1:length(eligible))]]
		rec = (r % 10 == 0) ? 0.5 : 10.0
		res = sirdif(alst, vlst, [seed_id], 0.33, WeibullDuration(rec, 2.0), 200, 0.5,
					-0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
					immunity=WeibullDuration(30.0, 2.0), rng=rng)
		daily = expand_daily_log(res["infection_log"], 200)
		t = Float64.(daily.time); C = daily.prop_cum_infected
		f = fit_richards(t, C; C0=C[1])
		inc = diff(C); ma7 = [sum(inc[k-3:k+3]) / 7 for k in 4:length(inc)-3]
		t_start = f.fit_ok ? richards_start(f.r_growth, f.t_i, f.a) : NaN
		push!(rows, (run=r, status=f.status, K=f.K, r=f.r_growth, t_i=f.t_i, a=f.a,
					frac=f.frac_at_inflection, t_start=t_start, fit_rmse=f.fit_rmse,
					obs_inc_peak=argmax(inc) - 0.5, obs_inc_peak_ma7=argmax(ma7) + 2.5,
					obs_peak_day=argmax(daily.n_infected) - 1,
					extinguished=any(daily.n_infected[2:end] .== 0),
					late_prev=mean(daily.n_infected[151:201]) / n))
	end
	show(DataFrame(rows); allrows=true)

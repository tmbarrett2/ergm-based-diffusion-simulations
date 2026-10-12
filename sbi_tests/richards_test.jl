#Richards Module Tests: UT-1.1, UT-1.2, and UT-1.3
#Jonathan H. Morgan, Ph.D.
#12 October 2026
#Revised 12 October 2026: t_i is the inflection time (the day of peak incidence);
#peak_incidence replaces richards_peak; fit_rmse is tested. The shape a now runs
#from -1 to 10: negative, Gompertz (a = 0), and monomolecular (a = -1) curves are
#tested, along with richards_start and inflection_fraction.

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

#   Checking for LsqFit (added to julia_env for richards.jl)
	if Base.find_package("LsqFit") === nothing
		error("LsqFit is not in this environment. Run: using Pkg; Pkg.add(\"LsqFit\")")
	end

################
#   PACKAGES   #
################
using DataFrames
using Random
using Statistics

#   Richards Module
	include("/workspace/ergm-based-diffusion-simulations/julia_env/src/richards.jl")
	using .richards_model

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Richards Tests: report a dictionary of checks
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
				println("  $(rpad(name, 34)) $(passed ? "PASS" : "FAIL")")
			end
			all_passed = all(values(checks))
			println("RESULT: $(all_passed ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")")

		#	Assembling Result
			return all_passed
	end

#	Helper Function for the Richards Tests: check that an expression throws a given error type
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

#	Helper Function for the Richards Tests: rise time from 10% to 90% of K
	function rise_width(r::Float64, a::Float64)
		"""
		Args:
			r::Float64: growth rate
			a::Float64: asymmetry exponent
		Returns:
			Float64: days for the Richards curve to rise from 0.1 K to 0.9 K
		Notes:
			Solving K [1 + a e^{-r(t - t_i)}]^{-1/a} = q K gives
			t_q = t_i - ln((q^{-a} - 1) / a) / r. The ratio is positive for
			a < 0 as well, so each log is taken of the ratio. For a = 0
			(Gompertz), t_q = t_i - ln(-ln q) / r.
		"""

		#	Gompertz case
			if a == 0.0
				return (log(-log(0.1)) - log(-log(0.9))) / r
			end

		#	Difference of the two crossing times
			return (log((0.1^(-a) - 1.0) / a) - log((0.9^(-a) - 1.0) / a)) / r
	end

#	Helper Function for the Richards Tests: prevalence series from a cumulative curve
	function prevalence_series(C::Vector{Float64}, D::Int)
		"""
		Args:
			C::Vector{Float64}: cumulative proportion on days 0, 1, ..., T
			D::Int: fixed infectious duration in days
		Returns:
			Vector{Float64}: proportion currently infected, I(t) = C(t) - C(t - D)
		Notes:
			Everyone infected by day t - D has recovered by day t. Before day 0
			the cumulative proportion is taken as 0.
		"""

		#	Difference of the cumulative curve over one infectious period
			n = length(C)
			I = Vector{Float64}(undef, n)
			for k in 1:n
				I[k] = C[k] - (k - D >= 1 ? C[k - D] : 0.0)
			end

		#	Assembling Result
			return I
	end

#	Helper Function for test_ut13_target: evaluate an extinct run under one policy
	function evaluate_extinct(region::TargetRegion, infected_t::Vector{Float64}, C0::Float64,
							  policy::Symbol, final_cum::Union{Float64, Nothing})
		"""
		Args:
			region::TargetRegion: acceptable outcome region
			infected_t::Vector{Float64}: daily proportion infected for the extinct run
			C0::Float64: known initial offset
			policy::Symbol: extinct_policy to apply
			final_cum::Union{Float64, Nothing}: final cumulative proportion, or nothing
		Returns:
			NamedTuple: the evaluate_target result
		Notes:
			An extinct run has no fitted curve, so the Richards parameters are NaN.
		"""

		#	Evaluate with the extinct status
			return evaluate_target(region, infected_t, NaN, NaN, NaN, NaN; C0=C0,
								   status=:extinct, extinct_policy=policy, final_cum=final_cum)
	end

#	simulate_richards: Offset, Seeds, and Validation
	function test_simulate_richards()
		"""
		Returns:
			Bool: true if every check passed
		Notes:
			Confirms the C0 offset, that seed = s draws the same numbers as
			rng = Xoshiro(s), that the global generator is left untouched,
			and that invalid inputs throw DomainError.
		"""

		#	Common inputs
			t = collect(0.0:1.0:200.0)
			checks = Dict{String, Bool}()

		#	Deterministic curve and offset
			clean = simulate_richards(t, 0.40, 0.15, 70.0, 1.0)
			shifted = simulate_richards(t, 0.40, 0.15, 70.0, 1.0; C0=0.02)
			checks["deterministic C equals C_true"] = clean.C == clean.C_true
			checks["offset adds C0 everywhere"] = all(isapprox.(shifted.C_true .- clean.C_true, 0.02; atol=1e-15))
			checks["params record C0"] = shifted.params.C0 == 0.02

		#	seed and rng give the same noise
			by_seed = simulate_richards(t, 0.40, 0.15, 70.0, 1.0; sigma=0.01, seed=42)
			by_rng = simulate_richards(t, 0.40, 0.15, 70.0, 1.0; sigma=0.01, rng=Xoshiro(42))
			checks["seed=42 matches rng=Xoshiro(42)"] = by_seed.C == by_rng.C

		#	Global generator is not reset by seed
			Random.seed!(7)
			x = rand()
			Random.seed!(7)
			simulate_richards(t, 0.40, 0.15, 70.0, 1.0; sigma=0.01, seed=42)
			y = rand()
			checks["global generator untouched"] = x == y

		#	Invalid inputs
			checks["K above 1 - C0 rejected"] = throws_error(() -> simulate_richards(t, 0.99, 0.15, 70.0, 1.0; C0=0.02), DomainError)
			checks["C0 = 1 rejected"] = throws_error(() -> simulate_richards(t, 0.40, 0.15, 70.0, 1.0; C0=1.0), DomainError)
			checks["a below -1 rejected"] = throws_error(() -> simulate_richards(t, 0.40, 0.15, 70.0, -1.5), DomainError)

		#	The extended shape range is accepted
			mono = simulate_richards(t, 0.40, 0.15, 70.0, -1.0)
			gomp = simulate_richards(t, 0.40, 0.15, 70.0, 0.0)
			checks["a = -1 and a = 0 accepted"] = all(isfinite, mono.C) && all(isfinite, gomp.C)

		#	Assembling Result
			return report_checks("simulate_richards: Offset, Seeds, and Validation", checks)
	end

#	Peak Timing, Peak Incidence, and Discretized Recovery Mean
	function test_peak_functions(; n_draws::Int=1_000_000, seed::Int=2026)
		"""
		Args:
			n_draws::Int: Monte Carlo draws per recovery distribution (default = 1_000_000)
			seed::Int: seed for the Monte Carlo generator (default = 2026)
		Returns:
			Bool: true if every check passed
		Notes:
			The inflection at t_i and peak_incidence are checked against a
			numerical derivative on a 0.001-day grid, for shapes from -1 to
			3. At a = -1 incidence is largest just after the start, t_i. The
			curve is checked against the classic form with location constant
			t_i + ln(a)/r (a > 0), the Gompertz curve (a = 0), and the
			monomolecular curve (a = -1), and for continuity through a = 0.
			richards_start and inflection_fraction are checked at known
			values. discrete_weibull_mean is checked against Monte
			Carlo draws made the way sirdif makes them:
			max(1, round(scale * E^(1/shape))) with E ~ Exp(1).
		"""

		#	Initialize checks
			checks = Dict{String, Bool}()

		#	Inflection and peak incidence against the numerical derivative
			tt = collect(0.0:0.001:250.0)
			for (K, r, ti, a) in [(0.40, 0.15, 70.0, 1.0), (0.60, 0.10, 90.0, 3.0), (0.25, 0.25, 50.0, 0.5),
								  (0.60, 0.09, 20.0, 1.0e-6), (0.60, 0.09, 20.0, 0.0),
								  (0.50, 0.10, 40.0, -0.5), (0.50, 0.10, 40.0, -1.0)]
				C = richards_curve(tt, [K, r, ti, a])
				dC = diff(C) ./ 0.001
				k = argmax(dC)
				checks["inflection at t_i, a=$(a)"] = abs(ti - (tt[k] + 0.0005)) < 0.01
				checks["peak incidence a=$(a)"] = abs(peak_incidence(K, r, a) - dC[k]) / dC[k] < 1e-4
			end

		#	Classic form with location constant t_i + ln(a)/r
			loc = 90.0 + log(3.0) / 0.1
			classic = 0.6 ./ (1.0 .+ exp.(-0.1 .* (tt .- loc))) .^ (1.0 / 3.0)
			checks["inflection form equals classic form"] = maximum(abs.(richards_curve(tt, [0.6, 0.1, 90.0, 3.0]) .- classic)) < 1e-12
			checks["logistic peak equals K r / 4"] = isapprox(peak_incidence(0.4, 0.15, 1.0), 0.4 * 0.15 / 4; rtol=1e-12)
			checks["Gompertz peak tends to K r / e"] = isapprox(peak_incidence(0.4, 0.15, 1.0e-6), 0.4 * 0.15 / exp(1.0); rtol=1e-5)
			checks["Gompertz peak equals K r / e"] = peak_incidence(0.4, 0.15, 0.0) == 0.4 * 0.15 * exp(-1.0)
			checks["monomolecular peak equals K r"] = peak_incidence(0.4, 0.15, -1.0) == 0.4 * 0.15

		#	Named members of the family and continuity through a = 0
			gompertz = 0.6 .* exp.(-exp.(-0.09 .* (tt .- 20.0)))
			checks["a = 0 equals the Gompertz curve"] = maximum(abs.(richards_curve(tt, [0.6, 0.09, 20.0, 0.0]) .- gompertz)) < 1e-14
			mono = 0.5 .* max.(0.0, 1.0 .- exp.(-0.1 .* (tt .- 40.0)))
			checks["a = -1 equals the monomolecular curve"] = maximum(abs.(richards_curve(tt, [0.5, 0.1, 40.0, -1.0]) .- mono)) < 1e-12
			near_zero = max(maximum(abs.(richards_curve(tt, [0.6, 0.09, 20.0, 1.0e-9]) .- gompertz)),
							maximum(abs.(richards_curve(tt, [0.6, 0.09, 20.0, -1.0e-9]) .- gompertz)))
			checks["continuous through a = 0"] = near_zero < 1e-7

		#	Start time of a curve with a < 0
			t0 = richards_start(0.1, 40.0, -0.5)
			checks["start time t_i + ln(-a)/r"] = isapprox(t0, 40.0 + log(0.5) / 0.1; rtol=1e-14)
			checks["curve is 0 before the start"] = richards_curve(t0 - 0.01, [0.5, 0.1, 40.0, -0.5]) == 0.0 &&
													richards_curve(t0 + 0.01, [0.5, 0.1, 40.0, -0.5]) > 0.0
			checks["no start for a >= 0"] = richards_start(0.1, 40.0, 0.0) == -Inf && richards_start(0.1, 40.0, 1.0) == -Inf

		#	Inflection fraction at known values and across the range
			checks["inflection fraction at a = 1, 0, -1"] = inflection_fraction(1.0) == 0.5 &&
															 inflection_fraction(0.0) == exp(-1.0) &&
															 inflection_fraction(-1.0) == 0.0
			checks["inflection fraction at a = -0.5"] = isapprox(inflection_fraction(-0.5), 0.25; rtol=1e-14)
			checks["inflection fraction increases with a"] = issorted([inflection_fraction(a) for a in -1.0:0.05:10.0])
			C_half = richards_curve(40.0, [0.5, 0.1, 40.0, -0.5])
			checks["inflection fraction equals y(t_i)/K"] = isapprox(C_half / 0.5, inflection_fraction(-0.5); rtol=1e-14)

		#	discrete_weibull_mean against Monte Carlo
			rng = Xoshiro(seed)
			for (lam, k) in [(10.0, 2.0), (3.0, 0.8), (14.0, 4.0)]
				draws = [max(1, round(Int, min(lam * randexp(rng)^(1.0 / k), 1.0e9))) for _ in 1:n_draws]
				checks["E[D] scale $(lam) shape $(k)"] = abs(discrete_weibull_mean(lam, k) - mean(draws)) < 0.02
			end
			checks["E[D] equals 1 at scale 0.5"] = abs(discrete_weibull_mean(0.5, 2.0) - 1.0) < 1e-3

		#	Termination for an extremely heavy tail
			value = 0.0
			elapsed = @elapsed (value = discrete_weibull_mean(10.0, 0.1))
			checks["heavy tail terminates"] = isfinite(value) && elapsed < 10.0

		#	peak_prevalence is peak incidence times E[D]
			expected = peak_incidence(0.4, 0.15, 1.3) * discrete_weibull_mean(10.0, 2.0)
			checks["peak_prevalence composition"] = peak_prevalence(0.4, 0.15, 1.3, 10.0, 2.0) == expected
			expected_neg = peak_incidence(0.4, 0.15, -0.5) * discrete_weibull_mean(10.0, 2.0)
			checks["peak_prevalence composition a < 0"] = peak_prevalence(0.4, 0.15, -0.5, 10.0, 2.0) == expected_neg

		#	Assembling Result
			return report_checks("Peak Timing, Peak Incidence, and Recovery Mean", checks)
	end

#	UT-1.1a: Noise-Free Recovery
	function test_ut11_deterministic(; verbose::Bool=true)
		"""
		Args:
			verbose::Bool: print the per-curve table (default = true)
		Returns:
			NamedTuple: (results::DataFrame, all_passed::Bool)
		Notes:
			The original 3 x 3 grid (t_i = 70, a = 1) plus the five UT-1.1
			parameter sets and four shapes at or below 0: a front-loaded
			curve (a = -0.5), a strongly front-loaded curve (a = -0.8), the
			Gompertz curve (a = 0), and the monomolecular curve (a = -1, the
			lower limit). Each is fitted with C0 = 0 and C0 = 0.02. On
			noise-free curves every parameter must be recovered to within
			1e-4 (relative; absolute for a = 0) and t_i to within 0.001 days,
			with status :ok.
		"""

		#	Build the cases
			cases = [(K, r, 70.0, 1.0) for K in [0.10, 0.40, 0.80] for r in [0.05, 0.15, 0.40]]
			append!(cases, [(0.40, 0.15, 70.0, 1.0), (0.80, 0.40, 30.0, 1.5), (0.60, 0.10, 90.0, 3.0),
							(0.25, 0.25, 50.0, 0.5), (0.10, 0.05, 90.0, 0.7),
							(0.60, 0.10, 30.0, -0.5), (0.30, 0.20, 25.0, -0.8),
							(0.60, 0.09, 20.0, 0.0), (0.50, 0.15, 2.0, -1.0)])
			t = collect(0.0:1.0:200.0)

		#	Pre-allocate results
			results = DataFrame(K = Float64[], r = Float64[], t_i = Float64[], a = Float64[], C0 = Float64[],
								status = Symbol[], K_err = Float64[], r_err = Float64[], ti_err = Float64[],
								a_err = Float64[], passed = Bool[])

		#	Fit each case with and without an offset
			for (K, r, ti, a) in cases, C0 in [0.0, 0.02]
				sim = simulate_richards(t, K, r, ti, a; C0=C0)
				fit = fit_richards(t, sim.C; C0=C0)
				K_err = abs(fit.K - K) / K
				r_err = abs(fit.r_growth - r) / r
				ti_err = abs(fit.t_i - ti)
				a_err = a == 0.0 ? abs(fit.a) : abs(fit.a - a) / abs(a)
				ok = fit.status == :ok && K_err < 1e-4 && r_err < 1e-4 && ti_err < 1e-3 && a_err < 1e-4
				push!(results, (K, r, ti, a, C0, fit.status, K_err, r_err, ti_err, a_err, ok))
			end

		#	Report
			all_passed = all(results.passed)
			if verbose
				println("=" ^ 60)
				println("UT-1.1a: Noise-Free Recovery ($(nrow(results)) curves)")
				println("=" ^ 60)
				println("  Passed: $(sum(results.passed)) / $(nrow(results))")
				if !all_passed
					println(results[.!results.passed, :])
				end
				println("RESULT: $(all_passed ? "ALL CURVES RECOVERED" : "SOME CURVES FAILED")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed)
	end

#	UT-1.1b: Noisy Recovery
	function test_ut11_noisy(; n_reps::Int=100, sigma::Float64=0.01, min_pass_rate::Float64=0.90,
							   verbose::Bool=true)
		"""
		Args:
			n_reps::Int: noisy replicates per case (default = 100)
			sigma::Float64: noise standard deviation (default = 0.01)
			min_pass_rate::Float64: share of replicates that must pass (default = 0.90)
			verbose::Bool: print the summary table (default = true)
		Returns:
			NamedTuple: (results::DataFrame, all_passed::Bool)
		Notes:
			Decision criterion (per replicate, fit_ok required): K within 5%,
			inflection time t_i within max(2 days, 5% of the 10-90% rise
			time), and peak incidence within 10%. The roadmap criterion on
			(K, r, t_i, a) is reported for information: under noise r and a
			trade off and are individually weakly identified.
			Five sets decide pass or fail, including a front-loaded curve
			(a = -0.5). The small, slow set (K = 0.10,
			r = 0.05) is reported for information: there sigma = 0.01 is 10%
			of the whole epidemic, an information limit rather than a
			fitting failure.
		"""

		#	Parameter sets: (K, r, t_i, a, decides)
			sets = [(0.40, 0.15, 70.0, 1.0, true),
					(0.80, 0.40, 30.0, 1.5, true),
					(0.60, 0.10, 90.0, 3.0, true),
					(0.25, 0.25, 50.0, 0.5, true),
					(0.60, 0.10, 30.0, -0.5, true),
					(0.10, 0.05, 90.0, 0.7, false)]
			t = collect(0.0:1.0:200.0)

		#	Pre-allocate results
			results = DataFrame(K = Float64[], r = Float64[], t_i = Float64[], a = Float64[], C0 = Float64[],
								decides = Bool[], fit_ok_rate = Float64[], derived_rate = Float64[],
								roadmap_rate = Float64[], statuses = String[], passed = Bool[])

		#	Iterate over sets and offsets
			case_id = 0
			for (K, r, ti, a, decides) in sets, C0 in [0.0, 0.02]
				#	True peak incidence and inflection tolerance
					case_id += 1
					pk_inc = peak_incidence(K, r, a)
					ti_tol = max(2.0, 0.05 * rise_width(r, a))

				#	Replicates
					n_ok = 0
					n_derived = 0
					n_roadmap = 0
					status_counts = Dict{Symbol, Int}()
					for rep in 1:n_reps
						sim = simulate_richards(t, K, r, ti, a; C0=C0, sigma=sigma, seed=10_000 * case_id + rep)
						fit = fit_richards(t, sim.C; C0=C0)
						status_counts[fit.status] = get(status_counts, fit.status, 0) + 1
						if fit.fit_ok
							n_ok += 1
							derived = abs(fit.K - K) / K < 0.05 &&
									  abs(fit.t_i - ti) < ti_tol &&
									  abs(fit.peak_incidence - pk_inc) / pk_inc < 0.10
							roadmap = abs(fit.K - K) / K < 0.05 && abs(fit.r_growth - r) / r < 0.10 &&
									  abs(fit.t_i - ti) < 2.0 && abs(fit.a - a) / abs(a) < 0.10
							n_derived += derived
							n_roadmap += roadmap
						end
					end

				#	Record the case
					status_str = join(["$(k)=$(v)" for (k, v) in sort(collect(status_counts))], ", ")
					derived_rate = n_derived / n_reps
					push!(results, (K, r, ti, a, C0, decides, n_ok / n_reps, derived_rate,
									n_roadmap / n_reps, status_str, derived_rate >= min_pass_rate))
			end

		#	Report
			all_passed = all(results.passed[results.decides])
			if verbose
				println("=" ^ 60)
				println("UT-1.1b: Noisy Recovery (sigma = $(sigma), $(n_reps) replicates per case)")
				println("=" ^ 60)
				show(results[:, [:K, :r, :t_i, :a, :C0, :decides, :fit_ok_rate, :derived_rate, :roadmap_rate, :passed]];
					 allrows=true, allcols=true)
				println()
				println("  Status counts by case:")
				for row in eachrow(results)
					println("    K=$(row.K) r=$(row.r) C0=$(row.C0): $(row.statuses)")
				end
				println("  Decision: derived_rate ≥ $(min_pass_rate) for the five deciding sets")
				println("RESULT: $(all_passed ? "ALL DECIDING CASES PASSED" : "SOME DECIDING CASES FAILED")")
			end

		#	Assembling Result
			return (results = results, all_passed = all_passed)
	end

#	UT-1.2: Degenerate Trajectories
	function test_ut12_degenerate()
		"""
		Returns:
			Bool: true if every check passed
		Notes:
			Every call is wrapped so that an exception escaping
			fit_richards counts as a failure of the no-exception check.
		"""

		#	Common inputs
			t = collect(0.0:1.0:200.0)
			checks = Dict{String, Bool}()
			escaped = false

		#	Fit every input, recording any exception that escapes
			results = Dict{String, Any}()
			inputs = Dict(
				"a: flat curve" => (t, fill(0.005, length(t)), 0.0),
				"b: small growth beyond C0" => (t, 0.02 .+ 0.015 ./ (1.0 .+ exp.(-0.1 .* (t .- 60.0))), 0.02),
				"b: same curve, C0 ignored" => (t, 0.02 .+ 0.015 ./ (1.0 .+ exp.(-0.1 .* (t .- 60.0))), 0.0),
				"c: non-starter" => (t, fill(1.0 / 614.0, length(t)), 1.0 / 614.0),
				"d: three points" => ([0.0, 1.0, 2.0], [0.0, 0.1, 0.2], 0.0),
				"e: vertical step" => (t, [ti < 100.0 ? 0.0 : 0.6 for ti in t], 0.0),
				"f: NaN in data" => (t, [k == 50 ? NaN : 0.4 / (1.0 + exp(-0.15 * (ti - 70.0))) for (k, ti) in enumerate(t)], 0.0),
				"g: mismatched lengths" => (t, fill(0.3, length(t) - 1), 0.0),
				"h: C0 outside [0, 1)" => (t, fill(0.3, length(t)), 1.0))
			for (name, (tt, C, C0)) in inputs
				try
					results[name] = fit_richards(tt, C; C0=C0)
				catch e
					escaped = true
					results[name] = nothing
				end
			end

		#	Expected classifications
			r = results
			checks["a: flat curve is extinct"] = r["a: flat curve"] !== nothing &&
												 r["a: flat curve"].extinct && !r["a: flat curve"].fit_ok &&
												 r["a: flat curve"].status == :extinct
			checks["b: growth below 0.02 beyond C0 is extinct"] = r["b: small growth beyond C0"] !== nothing &&
																  r["b: small growth beyond C0"].status == :extinct
			checks["b: without C0 it is not extinct"] = r["b: same curve, C0 ignored"] !== nothing &&
														r["b: same curve, C0 ignored"].status != :extinct
			checks["c: non-starter flagged, not viable"] = r["c: non-starter"] !== nothing &&
														   r["c: non-starter"].status == :nonstarter &&
														   !r["c: non-starter"].viable && !r["c: non-starter"].extinct
			checks["d: too few points"] = r["d: three points"] !== nothing &&
										  r["d: three points"].status == :too_few_points
			checks["e: vertical step rejected"] = r["e: vertical step"] !== nothing && !r["e: vertical step"].fit_ok &&
												  r["e: vertical step"].status in (:boundary, :too_steep, :no_convergence, :error)
			checks["f: NaN is invalid data"] = r["f: NaN in data"] !== nothing && r["f: NaN in data"].status == :invalid_data
			checks["g: mismatch is invalid data"] = r["g: mismatched lengths"] !== nothing &&
													r["g: mismatched lengths"].status == :invalid_data
			checks["h: bad C0 is invalid data"] = r["h: C0 outside [0, 1)"] !== nothing &&
												  r["h: C0 outside [0, 1)"].status == :invalid_data
			checks["no exception escapes fit_richards"] = !escaped

		#	Report the step-function status for the record
			if r["e: vertical step"] !== nothing
				println("  (vertical step status: $(r["e: vertical step"].status))")
			end

		#	Assembling Result
			return report_checks("UT-1.2: Degenerate Trajectories", checks)
	end

#	UT-1.3: Target Region Evaluation
	function test_ut13_target()
		"""
		Returns:
			Bool: true if every check passed
		Notes:
			The trajectory comes from a known Richards curve (K = 0.4,
			r = 0.15, t_i = 70, a = 1.3, C0 = 0.01) with a fixed 10-day
			infectious period, so I(t) = C(t) - C(t - 10). Bounds are set at
			1.1 times each actual quantity for the compliant region and at
			0.9 times one quantity for each single violation.
		"""

		#	Build the trajectory and its actual quantities
			t = collect(0.0:1.0:200.0)
			K, r, ti, a, C0 = 0.40, 0.15, 70.0, 1.3, 0.01
			C = simulate_richards(t, K, r, ti, a; C0=C0).C
			I = prevalence_series(C, 10)
			rec = (rec_scale=10.0, rec_shape=2.0)
			actual = (peak=maximum(I), surge=I[22], slope=maximum(diff(I)),
					  prev=peak_prevalence(K, r, a, rec.rec_scale, rec.rec_shape), final=C0 + K)
			checks = Dict{String, Bool}()

		#	Compliant region
			ok_region = TargetRegion(I_peak_max=1.1 * actual.peak, I_surge_max=1.1 * actual.surge,
									 S_max=1.1 * actual.slope, P_max_thresh=1.1 * actual.prev,
									 K_max=1.1 * actual.final)
			res = evaluate_target(ok_region, I, K, r, ti, a; C0=C0, rec...)
			checks["a: compliant run is in target"] = res.in_target === true && res.class == :established

		#	One violation at a time
			fields = [(:I_peak_max, :peak), (:I_surge_max, :surge), (:S_max, :slope),
					  (:P_max_thresh, :prevalence), (:K_max, :final_size)]
			bounds = Dict(:I_peak_max => actual.peak, :I_surge_max => actual.surge, :S_max => actual.slope,
						  :P_max_thresh => actual.prev, :K_max => actual.final)
			for (letter, (field, result_key)) in zip(["b", "c", "d", "e", "f"], fields)
				kw = Dict(f => 1.1 * v for (f, v) in bounds)
				kw[field] = 0.9 * bounds[field]
				region = TargetRegion(; kw...)
				res = evaluate_target(region, I, K, r, ti, a; C0=C0, rec...)
				others_hold = all(getfield(res, k) === true for (_, k) in fields if k != result_key)
				checks["$(letter): $(field) violation detected"] = res.in_target === false &&
																	getfield(res, result_key) === false && others_hold
			end

		#	Unconstrained region and configuration errors
			checks["g: TargetRegion() accepts any run"] = in_target(TargetRegion(), I, K, r, ti, a; C0=C0) === true
			checks["h: P_max without rec params throws"] = throws_error(() -> in_target(TargetRegion(P_max_thresh=0.05), I, K, r, ti, a), ArgumentError)
			checks["i: surge beyond series throws"] = throws_error(() -> in_target(TargetRegion(I_surge_max=0.1, surge_time=500), I, K, r, ti, a), ArgumentError)
			checks["j: unknown policy throws"] = throws_error(() -> in_target(TargetRegion(), I, K, r, ti, a; extinct_policy=:maybe), ArgumentError)
			checks["k: invalid bound rejected"] = throws_error(() -> TargetRegion(K_max=-0.1), DomainError)

		#	Non-starter: excluded under every policy
			for policy in (:exclude, :success, :failure)
				res = evaluate_target(ok_region, I, NaN, NaN, NaN, NaN; C0=C0, rec..., status=:nonstarter,
									  extinct_policy=policy)
				checks["l: non-starter excluded ($(policy))"] = res.in_target === nothing && res.class == :nonstarter
			end

		#	Extinct run under each policy
			I_ext = vcat([0.002, 0.003, 0.002, 0.001], zeros(197))
			lenient = TargetRegion(I_peak_max=0.10, K_max=0.20)
			checks["m: extinct :exclude gives nothing"] = evaluate_extinct(lenient, I_ext, C0, :exclude, 0.01).in_target === nothing
			checks["m: extinct :failure gives false"] = evaluate_extinct(lenient, I_ext, C0, :failure, 0.01).in_target === false
			checks["m: extinct :success within bounds"] = evaluate_extinct(lenient, I_ext, C0, :success, 0.01).in_target === true
			checks["m: extinct :success over K_max"] = evaluate_extinct(lenient, I_ext, C0, :success, 0.30).in_target === false
			checks["m: :success without final_cum throws"] = throws_error(() -> evaluate_extinct(lenient, I_ext, C0, :success, nothing), ArgumentError)

		#	Failed fits are excluded
			res_fail = evaluate_target(ok_region, I, 0.4, 5.0, 70.0, 1.0; C0=C0, rec..., status=:boundary)
			checks["n: failed fit excluded"] = res_fail.in_target === nothing && res_fail.class == :fit_failed
			res_nan = evaluate_target(ok_region, I, NaN, NaN, NaN, NaN; C0=C0, rec..., status=:ok)
			checks["o: :ok with NaN parameters excluded"] = res_nan.in_target === nothing && res_nan.class == :fit_failed

		#	Assembling Result
			return report_checks("UT-1.3: Target Region Evaluation", checks)
	end

#	Early-Peaking Epidemics: Regression for the Former Location Bound
	function test_early_peak_regression(; n_reps::Int=20, sigma::Float64=0.003)
		"""
		Args:
			n_reps::Int: noisy replicates per curve (default = 20)
			sigma::Float64: noise standard deviation (default = 0.003)
		Returns:
			Bool: true if every check passed
		Notes:
			Under the former parameterization, 75 of 100 simulated village
			epidemics failed with the location constant on its lower bound
			of 0. These curves use fitted values from three of those runs
			(location constant t_i + ln(a)/r at or below 0). Each must fit
			with status :ok, K within 2%, and t_i within 1 day, noise-free
			and in at least 90% of noisy replicates.
			Two exact curves check the former lower bound on a (1e-6): a
			Gompertz curve, which must now fit with a near 0 as an interior
			value, and a monomolecular curve starting on day 2, which must
			fit with a at its natural limit of -1.
		"""

		#	Curves: (K, r, t_i, a) from runs 1, 24, and 19 of the diagnostic
			t = collect(0.0:1.0:200.0)
			C0 = 1.0 / 614.0
			cases = [(0.627, 0.0874, 17.0, 0.226), (0.637, 0.0959, 20.1, 0.145), (0.648, 0.0865, 24.4, 0.121)]
			checks = Dict{String, Bool}()

		#	Richards cases, noise-free and noisy
			for (k, (K, r, ti, a)) in enumerate(cases)
				loc = ti + log(a) / r
				clean = simulate_richards(t, K, r, ti, a; C0=C0)
				fit = fit_richards(t, clean.C; C0=C0)
				checks["case $(k) noise-free (location $(round(loc, digits=2)))"] =
					fit.status == :ok && abs(fit.K - K) / K < 0.02 && abs(fit.t_i - ti) < 1.0
				n_ok = 0
				for rep in 1:n_reps
					noisy = simulate_richards(t, K, r, ti, a; C0=C0, sigma=sigma, seed=50_000 + 100 * k + rep)
					f = fit_richards(t, noisy.C; C0=C0)
					n_ok += f.status == :ok && abs(f.K - K) / K < 0.02 && abs(f.t_i - ti) < 1.0
				end
				checks["case $(k) noisy ($(n_ok)/$(n_reps))"] = n_ok >= 0.9 * n_reps
			end

		#	Exact Gompertz curve: K exp(-exp(-r (t - t_i)))
			Kg, rg, tig = 0.60, 0.09, 20.0
			gomp = C0 .+ Kg .* exp.(-exp.(-rg .* (t .- tig)))
			fg = fit_richards(t, gomp; C0=C0)
			checks["Gompertz curve fits across a = 0 (a = $(round(fg.a, sigdigits=3)))"] =
				fg.status == :ok && abs(fg.a) < 0.01 && abs(fg.K - Kg) / Kg < 0.02 && abs(fg.t_i - tig) < 1.0

		#	Exact monomolecular curve: K max(0, 1 - exp(-r (t - t_i)))
			Km, rm, tim = 0.50, 0.15, 2.0
			mono = C0 .+ Km .* max.(0.0, 1.0 .- exp.(-rm .* (t .- tim)))
			fm = fit_richards(t, mono; C0=C0)
			checks["monomolecular curve fits at a = -1 (a = $(round(fm.a, sigdigits=3)))"] =
				fm.status == :ok && abs(fm.a + 1.0) < 0.01 && abs(fm.K - Km) / Km < 0.02 && abs(fm.t_i - tim) < 1.0

		#	Assembling Result
			return report_checks("Early-Peaking Epidemics: Regression for the Former Location Bound", checks)
	end

#	fit_rmse: Residual Error of the Fit
	function test_fit_rmse(; n_reps::Int=50, sigma::Float64=0.01)
		"""
		Args:
			n_reps::Int: noisy replicates (default = 50)
			sigma::Float64: noise standard deviation (default = 0.01)
		Returns:
			Bool: true if every check passed
		Notes:
			A noise-free fit has a residual near 0. With C0 = 0.02 the curve
			stays two standard deviations above 0, so clamping at 0 is rare
			and the median residual should lie near sigma (0.8 to 1.2 times
			sigma). Runs with no completed fit report NaN.
		"""

		#	Common inputs
			t = collect(0.0:1.0:200.0)
			K, r, ti, a, C0 = 0.40, 0.15, 70.0, 1.3, 0.02
			checks = Dict{String, Bool}()

		#	Noise-free fit
			clean = fit_richards(t, simulate_richards(t, K, r, ti, a; C0=C0).C; C0=C0)
			checks["noise-free residual below 1e-5"] = clean.fit_ok && clean.fit_rmse < 1e-5

		#	Noisy fits
			rmse = Float64[]
			for rep in 1:n_reps
				f = fit_richards(t, simulate_richards(t, K, r, ti, a; C0=C0, sigma=sigma, seed=70_000 + rep).C; C0=C0)
				f.fit_ok && push!(rmse, f.fit_rmse)
			end
			med = isempty(rmse) ? NaN : median(rmse)
			checks["noisy median residual near sigma ($(round(med, digits=4)))"] = 0.8 * sigma <= med <= 1.2 * sigma

		#	No completed fit
			ext = fit_richards(t, fill(0.005, length(t)); C0=0.0)
			ns = fit_richards(t, fill(0.02, length(t)); C0=0.02)
			checks["extinct and non-starter report NaN"] = isnan(ext.fit_rmse) && isnan(ns.fit_rmse)

		#	Assembling Result
			return report_checks("fit_rmse: Residual Error of the Fit", checks)
	end

#############
#   TESTS   #
#############

#	Run Each Test
	simulate_passed = test_simulate_richards()
	peak_passed = test_peak_functions()
	ut11a = test_ut11_deterministic()
	ut11b = test_ut11_noisy()
	ut12_passed = test_ut12_degenerate()
	ut13_passed = test_ut13_target()
	early_passed = test_early_peak_regression()
	rmse_passed = test_fit_rmse()

#	Summary
	println("=" ^ 60)
	println("Richards Test Summary")
	println("=" ^ 60)
	println("  simulate_richards:          $(simulate_passed ? "PASS" : "FAIL")")
	println("  Peak and recovery functions: $(peak_passed ? "PASS" : "FAIL")")
	println("  UT-1.1a noise-free recovery: $(ut11a.all_passed ? "PASS" : "FAIL")")
	println("  UT-1.1b noisy recovery:      $(ut11b.all_passed ? "PASS" : "FAIL")")
	println("  UT-1.2 degenerate cases:     $(ut12_passed ? "PASS" : "FAIL")")
	println("  UT-1.3 target region:        $(ut13_passed ? "PASS" : "FAIL")")
	println("  Early-peak regression:       $(early_passed ? "PASS" : "FAIL")")
	println("  fit_rmse:                    $(rmse_passed ? "PASS" : "FAIL")")


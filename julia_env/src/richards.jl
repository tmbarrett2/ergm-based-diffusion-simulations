#Richards Model Simulation Extension
#Jonathan H. Morgan, Ph.D.
#12 April 2026
#Revised 12 October 2026: known initial offset C(0) (Option A), convergence,
#boundary, and transition checks, fit status codes, peak timing and incidence,
#and an explicit random number generator for simulate_richards. Target region,
#analytic peak prevalence, and in_target implemented (UT-1.3).

module richards_model

using LsqFit
using Random

# ============================================================
#	SIMULATION
# ============================================================

#	Simulate Richards Growth Curve
	function simulate_richards(t::Vector{Float64}, K::Float64, r::Float64, t_i::Float64, a::Float64;
							   C0::Float64=0.0,
							   sigma::Union{Float64, Nothing}=nothing, seed::Union{Int, Nothing}=nothing,
							   rng::Union{AbstractRNG, Nothing}=nothing)
		"""
		Args:
			t::Vector{Float64}: time points at which to evaluate the curve
			K::Float64: final epidemic size beyond the initial offset (0 < K ≤ 1 - C0)
			r::Float64: per-capita growth rate (r > 0)
			t_i::Float64: inflection time parameter
			a::Float64: asymmetry exponent (a > 0; a = 1 recovers logistic)
			C0::Float64: known initial offset, the seed fraction (default = 0.0)
			sigma::Union{Float64, Nothing}: standard deviation of additive Gaussian noise (default = nothing, deterministic)
			seed::Union{Int, Nothing}: seed for a local Xoshiro generator when sigma is set (default = nothing)
			rng::Union{AbstractRNG, Nothing}: generator for the noise when seed is not set (default = nothing, global generator)
		Returns:
			NamedTuple with fields:
				t::Vector{Float64}: input time vector
				C::Vector{Float64}: cumulative incidence values (noisy if sigma is set)
				C_true::Vector{Float64}: deterministic curve C0 + Richards (always clean)
				params::NamedTuple: (K, r, t_i, a, C0) used to generate the curve
		Notes:
			C(t) = C0 + K / [1 + exp(-r * (t - t_i))]^(1/a)
			C0 mirrors the simulator, which starts at seeds/N rather than 0.
			The seed keyword now builds a local Xoshiro(seed) instead of calling
			Random.seed!(seed). Xoshiro(s) reproduces Random.seed!(s) exactly
			(confirmed in UT-1.7), so earlier seeded curves are unchanged, but
			the global generator is no longer reset as a side effect.
		"""

		#	Validation
			if C0 < 0.0 || C0 >= 1.0
				throw(DomainError(C0, "C0 must be in [0, 1)"))
			end
			if K <= 0.0 || K > 1.0 - C0
				throw(DomainError(K, "K must be in (0, 1 - C0]"))
			end
			if r <= 0.0
				throw(DomainError(r, "r must be positive"))
			end
			if a <= 0.0
				throw(DomainError(a, "a must be positive"))
			end
			if length(t) < 2
				throw(ArgumentError("t must contain at least 2 time points"))
			end
			if sigma !== nothing && sigma < 0.0
				throw(DomainError(sigma, "sigma must be non-negative"))
			end

		#	Compute deterministic curve: offset plus Richards growth
			C_true = C0 .+ richards_curve(t, [K, r, t_i, a])

		#	Add optional Gaussian noise
			if sigma !== nothing && sigma > 0.0
				#	Choose the generator: seed first, then rng, then the global generator
					noise_rng = seed !== nothing ? Xoshiro(seed) :
								rng !== nothing ? rng : Random.default_rng()

				#	Generate noisy curve clamped to [0, 1]
					noise = sigma .* randn(noise_rng, length(t))
					C = clamp.(C_true .+ noise, 0.0, 1.0)
			else
				#	Deterministic case
					C = copy(C_true)
			end

		#	Assembling Result
			return (t=t, C=C, C_true=C_true, params=(K=K, r=r, t_i=t_i, a=a, C0=C0))
	end
	@doc raw"""
	**Description**
	Generate a Richards growth curve with known parameters and an optional initial offset, optionally perturbed by additive Gaussian noise. Used to produce synthetic cumulative incidence trajectories for testing `fit_richards` parameter recovery.

	**Usage**
	```julia
	simulate_richards(t, K, r, t_i, a; C0=0.0, sigma=nothing, seed=nothing, rng=nothing)
	```

	**Arguments**
	- `t::Vector{Float64}`: Time points at which to evaluate the curve.
	- `K::Float64`: Final epidemic size beyond the offset, as a proportion (0 < K ≤ 1 - C0).
	- `r::Float64`: Per-capita growth rate (must be positive).
	- `t_i::Float64`: Inflection time parameter.
	- `a::Float64`: Asymmetry exponent (must be positive). `a = 1` recovers the standard logistic curve.
	- `C0::Float64`: Known initial offset (default `0.0`). In simulator output this is the seed fraction, seeds/N.
	- `sigma::Union{Float64, Nothing}`: Standard deviation of additive Gaussian noise (default `nothing` for deterministic output).
	- `seed::Union{Int, Nothing}`: Seed for a local `Xoshiro` generator when noise is applied (default `nothing`).
	- `rng::Union{AbstractRNG, Nothing}`: Generator for the noise when `seed` is not set (default `nothing`, the global generator).

	**Details**
	The curve is:

	``C(t) = C_0 + \frac{K}{[1 + \exp(-r \cdot (t - t_i))]^{1/a}}``

	When `sigma` is `nothing` or `0.0`, the output is the exact deterministic curve. When `sigma > 0`, independent Gaussian noise ``N(0, \sigma^2)`` is added to each point and the result is clamped to `[0.0, 1.0]` since `C` represents a proportion.

	A `seed` builds a local `Xoshiro(seed)`, which draws the same numbers as `Random.seed!(seed)` without resetting the global generator.

	The deterministic curve is always returned in `C_true` regardless of the noise setting, enabling direct comparison between noisy observations and ground truth.

	**Value**
	A `NamedTuple` with fields:
	- `t::Vector{Float64}`: The input time vector.
	- `C::Vector{Float64}`: Cumulative incidence values (noisy if `sigma` is set, otherwise identical to `C_true`).
	- `C_true::Vector{Float64}`: The clean deterministic curve.
	- `params::NamedTuple`: The generating parameters `(K, r, t_i, a, C0)`.

	**Examples**
	```julia
	using Random

	#	Define time grid
		t = collect(0.0:1.0:200.0)

	#	Deterministic curve
		result = simulate_richards(t, 0.15, 0.10, 70.0, 1.3)
		println("Final size: ", result.C[end])

	#	Noisy curve with three seeds on a 150-person village
		noisy = simulate_richards(t, 0.15, 0.10, 70.0, 1.3; C0=3/150, sigma=0.01, seed=42)
		println("Max deviation: ", maximum(abs.(noisy.C .- noisy.C_true)))
	```

	**References**
	Richards, F.J. (1959). A flexible growth function for empirical use.
	*Journal of Experimental Botany*, 10(2), 290–301.
	""" simulate_richards

# ============================================================
#	FITTING
# ============================================================

#	Helper Function for fit_richards: Richards curve
	function richards_curve(t, p)
		"""
		Args:
			t: time points (vector or scalar)
			p: parameter vector [K, r, t_i, a]
		Returns:
			Vector{Float64} (or Float64): K / [1 + exp(-r (t - t_i))]^(1/a)
		Notes:
			Model function passed to LsqFit.curve_fit. Moved out of
			fit_richards to module scope; the formula is unchanged.
			Overflow of exp() yields Inf in the denominator and 0.0 in the
			result, which is the correct limit.
		"""

		#	Unpack and evaluate
			K, r, ti, a = p
			return K ./ (1.0 .+ exp.(-r .* (t .- ti))) .^ (1.0 / a)
	end

#	Helper Function for fit_richards: test whether a parameter sits on a bound
	function at_bound(x::Float64, lo::Float64, hi::Float64; rtol::Float64=1e-6)
		"""
		Args:
			x::Float64: fitted parameter value
			lo::Float64: lower bound
			hi::Float64: upper bound
			rtol::Float64: tolerance as a fraction of hi - lo (default = 1e-6)
		Returns:
			Bool: true if x is within rtol * (hi - lo) of either bound
		"""

		#	Compare with both bounds
			tol = rtol * (hi - lo)
			return (x - lo) <= tol || (hi - x) <= tol
	end

#	Richards Peak Timing and Incidence
	function richards_peak(K::Float64, r::Float64, t_i::Float64, a::Float64)
		"""
		Args:
			K::Float64: final size (beyond the offset)
			r::Float64: growth rate
			t_i::Float64: inflection time parameter
			a::Float64: asymmetry exponent
		Returns:
			NamedTuple: (t_peak::Float64, peak_incidence::Float64)
		Notes:
			dC/dt peaks at t_peak = t_i - ln(a)/r with value
			K * r * (1 + a)^(-(1 + a)/a). t_peak equals t_i only when a = 1.
			These corrected expressions replace the formula previously
			quoted in the roadmap (§2.3), which holds only for a = 1.
		"""

		#	Peak time and peak incidence
			t_peak = t_i - log(a) / r
			peak_incidence = K * r * (1.0 + a)^(-(1.0 + a) / a)

		#	Assembling Result
			return (t_peak=t_peak, peak_incidence=peak_incidence)
	end
	@doc raw"""
	**Description**
	Time and height of peak incidence for a Richards curve.

	**Usage**
	`richards_peak(K, r, t_i, a)`

	**Details**
	For ``C(t) = K[1 + e^{-r(t - t_i)}]^{-1/a}``, the incidence ``\dot C(t)`` is maximal at

	``t^* = t_i - \frac{\ln a}{r}, \qquad \dot C(t^*) = K r (1 + a)^{-(1 + a)/a}.``

	The parameter ``t_i`` is the peak time only when ``a = 1``, where ``\dot C(t^*) = Kr/4``. Under noise, ``r``, ``t_i``, and ``a`` are individually weakly identified, but ``t^*`` and ``\dot C(t^*)`` are recovered reliably, so they are reported with every fit.

	**Value**
	`NamedTuple` `(t_peak, peak_incidence)`.

	**See Also**
	`fit_richards`
	""" richards_peak

#	Helper Function for fit_richards: assemble a fit result
	function fit_result(K::Float64, r::Float64, ti::Float64, a::Float64,
						fit_ok::Bool, extinct::Bool, status::Symbol)
		"""
		Args:
			K, r, ti, a::Float64: fitted (or NaN) Richards parameters
			fit_ok::Bool: true only for a usable fit
			extinct::Bool: true if the run was classified as extinct
			status::Symbol: reason code (see fit_richards)
		Returns:
			NamedTuple: K, r_growth, t_i, a, t_peak, peak_incidence, fit_ok, extinct, status
		Notes:
			Replaces the nested result constructor that used to sit inside
			fit_richards. Peak timing and incidence are computed whenever
			all four parameters are finite.
		"""

		#	Derived peak quantities
			if all(isfinite, (K, r, ti, a)) && r > 0.0 && a > 0.0
				pk = richards_peak(K, r, ti, a)
				t_peak = pk.t_peak
				peak_inc = pk.peak_incidence
			else
				t_peak = NaN
				peak_inc = NaN
			end

		#	Assembling Result
			return (K=K, r_growth=r, t_i=ti, a=a, t_peak=t_peak, peak_incidence=peak_inc,
					fit_ok=fit_ok, extinct=extinct, status=status)
	end

#	Fit Richards Growth Model
	function fit_richards(t::Vector{Float64}, C::Vector{Float64};
						  C0::Float64=0.0,
						  extinction_threshold::Float64=0.02,
						  min_points::Int=5,
						  min_transition_points::Int=3,
						  bound_rtol::Float64=1e-6,
						  max_iter::Int=2000)
		"""
		Args:
			t::Vector{Float64}: time points
			C::Vector{Float64}: cumulative proportion ever infected
			C0::Float64: known initial offset, seeds/N for simulator output (default = 0.0)
			extinction_threshold::Float64: C[end] - C0 below this → extinct (default = 0.02)
			min_points::Int: minimum number of data points required (default = 5)
			min_transition_points::Int: minimum observations between 10% and 90% of the fitted K (default = 3)
			bound_rtol::Float64: tolerance for the boundary check, as a fraction of each bound range (default = 1e-6)
			max_iter::Int: iteration limit passed to LsqFit (default = 2000)
		Returns:
			NamedTuple with fields:
				K::Float64: fitted final size beyond C0
				r_growth::Float64: fitted per-capita growth rate
				t_i::Float64: fitted inflection time parameter
				a::Float64: fitted asymmetry exponent
				t_peak::Float64: time of peak incidence, t_i - ln(a)/r
				peak_incidence::Float64: peak of dC/dt, K r (1 + a)^(-(1 + a)/a)
				fit_ok::Bool: true only for a usable fit
				extinct::Bool: true if C[end] - C0 < extinction_threshold
				status::Symbol: :ok, :extinct, :too_few_points, :invalid_data,
								:no_convergence, :boundary, :too_steep,
								:implausible, or :error
		Notes:
			Option A (open decision 1, 12 October 2026): the known offset C0
			is subtracted before fitting, so K is the share infected beyond
			the seeds and extinction is judged on growth beyond the seeds.
			No exception propagates out of this function; every failure
			returns fit_ok = false with a status code.
		"""

		#	Check data length and values
			if length(t) != length(C) || length(t) < min_points
				status = length(t) != length(C) ? :invalid_data : :too_few_points
				return fit_result(NaN, NaN, NaN, NaN, false, false, status)
			end
			if !all(isfinite, t) || !all(isfinite, C) || !(0.0 <= C0 < 1.0)
				return fit_result(NaN, NaN, NaN, NaN, false, false, :invalid_data)
			end

		#	Subtract the known offset (Option A)
			y = C .- C0

		#	Check for extinction: growth beyond the seeds is below threshold
			if y[end] < extinction_threshold
				return fit_result(NaN, NaN, NaN, NaN, false, true, :extinct)
			end

		#	Compute initial parameter guesses from data
		#	t_i0 uses the first crossing of half the final size, which is
		#	less sensitive to noise than the largest single-day increase
			K0 = y[end]
			r0 = 0.3
			half_idx = findfirst(>=(K0 / 2.0), y)
			t_i0 = half_idx === nothing ? t[div(length(t), 2)] : t[half_idx]
			a0 = 1.0
			p0 = [K0, r0, t_i0, a0]

		#	Set parameter bounds; K cannot exceed the susceptible share 1 - C0
			T_max = t[end]
			lower = [0.0,      0.001, 0.0,   0.1]
			upper = [1.0 - C0, 5.0,   T_max, 10.0]

		#	Clamp initial guesses to bounds
			p0 = clamp.(p0, lower, upper)

		#	Attempt curve fitting
			try
				#	Run nonlinear least squares
					fit = curve_fit(richards_curve, t, y, p0; lower=lower, upper=upper, maxIter=max_iter)

				#	Extract fitted parameters
					K_fit = fit.param[1]
					r_fit = fit.param[2]
					ti_fit = fit.param[3]
					a_fit = fit.param[4]

				#	Check for non-finite parameters
					if !all(isfinite, fit.param)
						return fit_result(NaN, NaN, NaN, NaN, false, false, :error)
					end

				#	Check convergence: LsqFit returns rather than throws when it stops early
					if !fit.converged
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :no_convergence)
					end

				#	Check for implausible fit
					if K_fit < extinction_threshold
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :implausible)
					end

				#	Check bounds on r, t_i, and a; K may legitimately reach 1 - C0
					if at_bound(r_fit, lower[2], upper[2]; rtol=bound_rtol) ||
					   at_bound(ti_fit, lower[3], upper[3]; rtol=bound_rtol) ||
					   at_bound(a_fit, lower[4], upper[4]; rtol=bound_rtol)
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :boundary)
					end

				#	Check that the rise is resolved by the time grid
				#	A near-vertical step leaves no observations between 10% and 90% of K
					n_transition = count(v -> 0.1 * K_fit <= v <= 0.9 * K_fit, y)
					if n_transition < min_transition_points
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :too_steep)
					end

				#	Assembling Result
					return fit_result(K_fit, r_fit, ti_fit, a_fit, true, false, :ok)

			catch e
				#	Numerical failure inside LsqFit
					return fit_result(NaN, NaN, NaN, NaN, false, false, :error)
			end
	end
	@doc raw"""
	**Description**
	Fit the Richards growth model to cumulative incidence data using nonlinear least squares, after subtracting the known initial offset. Returns fitted parameters, peak timing and incidence, and diagnostic flags for degenerate trajectories.

	**Usage**
	```julia
	fit_richards(t, C; C0=0.0, extinction_threshold=0.02, min_points=5,
				 min_transition_points=3, bound_rtol=1e-6, max_iter=2000)
	```

	**Arguments**
	- `t::Vector{Float64}`: Time points of the observations.
	- `C::Vector{Float64}`: Cumulative proportion ever infected at each time point.
	- `C0::Float64`: Known initial offset (default `0.0`). For simulator output, pass the seed fraction, seeds/N.
	- `extinction_threshold::Float64`: If `C[end] - C0` falls below this value, the run is classified as extinct and fitting is skipped (default `0.02`).
	- `min_points::Int`: Minimum number of data points required for fitting (default `5`).
	- `min_transition_points::Int`: Minimum number of observations between 10% and 90% of the fitted `K` (default `3`).
	- `bound_rtol::Float64`: Tolerance for the boundary check, as a fraction of each parameter's bound range (default `1e-6`).
	- `max_iter::Int`: Iteration limit passed to `LsqFit` (default `2000`).

	**Details**
	The model fitted to ``y(t) = C(t) - C_0`` is:

	``y(t) = \frac{K}{[1 + \exp(-r \cdot (t - t_i))]^{1/a}}``

	Subtracting the known offset (Option A) keeps the fit at four parameters and makes `K` the share of the population infected beyond the seeds. The extinction rule then works for any number of seed nodes.

	Initial parameter guesses are derived from the data:
	- ``K_0 = y[\text{end}]``
	- ``r_0 = 0.3``
	- ``t_{i,0}`` = first time ``y`` reaches ``K_0 / 2``
	- ``a_0 = 1.0``

	Parameter bounds enforced during fitting:

	| Parameter | Lower | Upper | Notes |
	|-----------|-------|-------|-------|
	| K | 0.0 | 1 - C0 | Share infected beyond the seeds |
	| r | 0.001 | 5.0 | Growth rate; upper bound prevents divergence |
	| t_i | 0.0 | T_max | Inflection parameter within the horizon |
	| a | 0.1 | 10.0 | Shape parameter; a=1 is logistic |

	Degenerate trajectory handling (`status` in parentheses):
	- **Too few points** (`:too_few_points`): `length(t) < min_points`.
	- **Invalid data** (`:invalid_data`): mismatched lengths, non-finite values, or `C0` outside `[0, 1)`.
	- **Extinction** (`:extinct`): `C[end] - C0 < extinction_threshold`; `extinct=true`.
	- **Non-convergence** (`:no_convergence`): `LsqFit` stopped without meeting its convergence tolerances. `LsqFit` reports this through `converged = false` rather than by throwing.
	- **Implausible fit** (`:implausible`): fitted `K < extinction_threshold`.
	- **Boundary fit** (`:boundary`): `r`, `t_i`, or `a` on a bound. `K` may reach `1 - C0`, which means everyone beyond the seeds was infected.
	- **Too steep** (`:too_steep`): fewer than `min_transition_points` observations between 10% and 90% of the fitted `K`, so the rise is not resolved by the time grid.
	- **Numerical error** (`:error`): `LsqFit` threw, or returned non-finite parameters.

	No exception propagates out of `fit_richards`. When `fit_ok=false`, parameters are `NaN` for failures that occur before or during fitting, and the fitted values are returned for diagnostic inspection when a fit completed but failed a check.

	Under noise, `r`, `t_i`, and `a` trade off against one another and are individually weakly identified. `t_peak` and `peak_incidence` are recovered reliably and are reported with every fit (see `richards_peak`).

	**Value**
	A `NamedTuple` with fields:
	- `K::Float64`: Fitted final size beyond `C0`.
	- `r_growth::Float64`: Fitted per-capita growth rate.
	- `t_i::Float64`: Fitted inflection time parameter.
	- `a::Float64`: Fitted asymmetry exponent.
	- `t_peak::Float64`: Time of peak incidence.
	- `peak_incidence::Float64`: Peak of ``\dot C``.
	- `fit_ok::Bool`: `true` if fitting succeeded and passed every check.
	- `extinct::Bool`: `true` if the run was classified as extinct.
	- `status::Symbol`: Reason code, `:ok` for a usable fit.

	**Examples**
	```julia
	#	Generate a known curve with three seeds on 150 nodes
		t = collect(0.0:1.0:200.0)
		result = simulate_richards(t, 0.40, 0.15, 70.0, 1.0; C0=3/150)

	#	Recover parameters, passing the known offset
		fit = fit_richards(t, result.C; C0=3/150)
		println("K: ", fit.K, "  r: ", fit.r_growth, "  t_i: ", fit.t_i, "  a: ", fit.a)
		println("t_peak: ", fit.t_peak, "  status: ", fit.status)
	```

	**References**
	Richards, F.J. (1959). A flexible growth function for empirical use.
	*Journal of Experimental Botany*, 10(2), 290–301.

	**See Also**
	`simulate_richards`, `richards_peak`, `LsqFit.curve_fit`
	""" fit_richards

# ============================================================
#	TARGET REGION
# ============================================================

#	Target Region Struct
	struct TargetRegion
		I_peak_max   :: Union{Float64, Nothing}
		I_surge_max  :: Union{Float64, Nothing}
		surge_time   :: Int
		S_max        :: Union{Float64, Nothing}
		P_max_thresh :: Union{Float64, Nothing}
		K_max        :: Union{Float64, Nothing}
	end

#	Target Region Keyword Constructor
	function TargetRegion(; I_peak_max::Union{Real, Nothing}=nothing,
							I_surge_max::Union{Real, Nothing}=nothing,
							surge_time::Integer=21,
							S_max::Union{Real, Nothing}=nothing,
							P_max_thresh::Union{Real, Nothing}=nothing,
							K_max::Union{Real, Nothing}=nothing)
		"""
		Args:
			I_peak_max: bound on max(infected_t), proportion (default = nothing, unconstrained)
			I_surge_max: bound on infected_t at surge_time, proportion (default = nothing)
			surge_time::Integer: day of the early-infection check (default = 21)
			S_max: bound on the largest daily increase in infected_t (default = nothing)
			P_max_thresh: bound on analytic peak prevalence (default = nothing)
			K_max: bound on final size, C0 + K (default = nothing)
		Returns:
			TargetRegion
		Notes:
			TargetRegion() gives an all-unconstrained region with surge_time = 21,
			as the roadmap specifies. Bounds must be positive and finite;
			surge_time must be non-negative.
		"""

		#	Validation of the bounds
			for (name, v) in (("I_peak_max", I_peak_max), ("I_surge_max", I_surge_max),
							  ("S_max", S_max), ("P_max_thresh", P_max_thresh), ("K_max", K_max))
				if v !== nothing && (!isfinite(v) || v <= 0.0)
					throw(DomainError(v, "$(name) must be positive and finite, or nothing"))
				end
			end
			if surge_time < 0
				throw(DomainError(surge_time, "surge_time must be ≥ 0"))
			end

		#	Convert to stored types
			return TargetRegion(I_peak_max === nothing ? nothing : Float64(I_peak_max),
								I_surge_max === nothing ? nothing : Float64(I_surge_max),
								Int(surge_time),
								S_max === nothing ? nothing : Float64(S_max),
								P_max_thresh === nothing ? nothing : Float64(P_max_thresh),
								K_max === nothing ? nothing : Float64(K_max))
	end
	@doc raw"""
	**Description**
	The planner-defined acceptable outcome region ``\mathcal{A}``: a conjunction of optional upper bounds on five epidemiological quantities.

	**Usage**
	```julia
	TargetRegion()
	TargetRegion(; I_peak_max=nothing, I_surge_max=nothing, surge_time=21,
				   S_max=nothing, P_max_thresh=nothing, K_max=nothing)
	```

	**Arguments**
	- `I_peak_max`: Peak infection size, the largest proportion simultaneously infected ("does this overwhelm capacity?").
	- `I_surge_max`: Early-infection size, the proportion infected on day `surge_time` ("how fast do we need to mobilize?").
	- `surge_time`: Day of the early-infection check (default 21).
	- `S_max`: Slope, the largest one-day increase in the proportion infected ("how fast is it accelerating?").
	- `P_max_thresh`: Prevalence, the analytic peak prevalence from the Richards fit and the recovery distribution ("what is the load on the health system?").
	- `K_max`: Final size, the total proportion ever infected, `C0 + K` ("what is the total impact?").

	**Details**
	Every bound is optional; `nothing` leaves that quantity unconstrained. `TargetRegion()` is fully unconstrained with `surge_time = 21`. Bounds must be positive and finite.

	**Value**
	A `TargetRegion`.

	**See Also**
	`in_target`, `evaluate_target`, `peak_prevalence`
	""" TargetRegion

#	Helper Function for peak_prevalence: mean of the discretized Weibull duration
	function discrete_weibull_mean(scale::Float64, shape::Float64; tol::Float64=1e-14,
	                               max_days::Int=1_000_000)
		"""
		Args:
			scale::Float64: Weibull scale λ (> 0)
			shape::Float64: Weibull shape k (> 0)
			tol::Float64: stop once a term falls below tol past the scale (default = 1e-14)
			max_days::Int: upper limit on the summation, in days (default = 1_000_000)
		Returns:
			Float64: E[D] for D = max(1, round(X)), X ~ Weibull(λ, k)
		Notes:
			This is the duration sirdif actually uses (draw_duration in
			diffusion_sim.jl). Because D ≥ 1 is integer valued,
			E[D] = Σ_{d≥1} P(D ≥ d) = 1 + Σ_{d≥2} exp(-((d - 0.5)/λ)^k).
			It differs from λ Γ(1 + 1/k) mainly for short durations:
			at λ = 0.5, k = 2 the continuous mean is 0.44 but E[D] = 1.00.
			The sum stops at max_days, which guarantees termination. The
			truncation matters only for extremely heavy tails (shape below
			about 0.2), far outside the prior on k_rec; for those it returns
			the mean of min(D, max_days).
		"""

		#	Validation
			if !(scale > 0.0) || !(shape > 0.0) || !isfinite(scale) || !isfinite(shape)
				throw(DomainError((scale, shape), "scale and shape must be positive and finite"))
			end

		#	Sum the survival terms
			total = 1.0
			d = 2
			while true
				term = exp(-((d - 0.5) / scale)^shape)
				total += term
				if (term < tol && d > scale) || d >= max_days
					break
				end
				d += 1
			end

		#	Assembling Result
			return total
	end

#	Analytic Peak Prevalence
	function peak_prevalence(K::Float64, r::Float64, a::Float64,
							 rec_scale::Float64, rec_shape::Float64)::Float64
		"""
		Args:
			K::Float64: fitted final size beyond C0
			r::Float64: fitted growth rate
			a::Float64: fitted asymmetry exponent
			rec_scale::Float64: Weibull recovery scale
			rec_shape::Float64: Weibull recovery shape
		Returns:
			Float64: P_max = I_max * E[D]
		Notes:
			I_max = K r (1 + a)^(-(1 + a)/a) is the peak incidence of the
			Richards curve (richards_peak; t_i does not affect its height).
			E[D] is the mean of the discretized recovery duration
			(discrete_weibull_mean). P_max approximates peak prevalence by
			treating incidence as constant near the peak for one mean
			infectious period; it overstates the true peak when the
			incidence peak is narrow relative to E[D].
		"""

		#	Peak incidence times mean infectious duration
			I_max = richards_peak(K, r, 0.0, a).peak_incidence
			return I_max * discrete_weibull_mean(rec_scale, rec_shape)
	end
	@doc raw"""
	**Description**
	Analytic peak prevalence from the Richards parameters and the Weibull recovery distribution, for prospective use without a simulation.

	**Usage**
	`peak_prevalence(K, r, a, rec_scale, rec_shape)`

	**Details**
	``P_{\max} = I_{\max} \times E[D]`` with ``I_{\max} = K r (1 + a)^{-(1+a)/a}`` (the Richards peak incidence) and ``E[D] = 1 + \sum_{d \ge 2} \exp\{-((d - 1/2)/\lambda)^k\}``, the mean of the recovery duration as `sirdif` draws it (rounded to whole days with a floor of one day). For ``a = 1``, ``I_{\max} = Kr/4``.

	The product treats incidence as constant near its peak for one mean infectious period, so it is an approximation. It overstates the true peak prevalence when the incidence peak is narrow relative to ``E[D]``.

	**Value**
	`Float64` proportion of the population.

	**See Also**
	`richards_peak`, `discrete_weibull_mean`, `in_target`
	""" peak_prevalence

#	Helper Function for evaluate_target: test one optional upper bound
	function within_bound(value::Float64, bound::Union{Float64, Nothing})
		"""
		Args:
			value::Float64: observed or fitted quantity
			bound::Union{Float64, Nothing}: upper bound, or nothing for unconstrained
		Returns:
			Union{Bool, Nothing}: nothing if unconstrained, otherwise value ≤ bound
			(a non-finite value fails)
		"""

		#	Unconstrained or compared
			if bound === nothing
				return nothing
			end
			return isfinite(value) && value <= bound
	end

#	Target Region Constraint Evaluation
	function evaluate_target(region::TargetRegion,
							 infected_t::Vector{Float64},
							 K::Float64, r::Float64, ti::Float64, a::Float64;
							 C0::Float64=0.0,
							 rec_scale::Union{Float64, Nothing}=nothing,
							 rec_shape::Union{Float64, Nothing}=nothing)
		"""
		Args:
			region::TargetRegion: acceptable outcome region
			infected_t::Vector{Float64}: proportion currently infected on days 0, 1, ..., T
										 (one value per day; summarize.jl expands the simulator log)
			K, r, ti, a::Float64: fitted Richards parameters (K beyond C0)
			C0::Float64: known initial offset, seeds/N (default = 0.0)
			rec_scale, rec_shape: Weibull recovery parameters, required when P_max_thresh is set
		Returns:
			NamedTuple: (peak, surge, slope, prevalence, final_size, in_target)
			Each constraint is true, false, or nothing (unconstrained).
		Notes:
			Throws ArgumentError if P_max_thresh is set without rec_scale and
			rec_shape, if infected_t is empty, or if surge_time lies beyond
			the series. These are configuration errors, not data outcomes.
			Non-finite Richards parameters (a failed or extinct fit) make
			in_target false, as the roadmap specifies.
		"""

		#	Validation of the configuration
			if region.P_max_thresh !== nothing && (rec_scale === nothing || rec_shape === nothing)
				throw(ArgumentError("P_max_thresh is set, so rec_scale and rec_shape must be provided"))
			end
			if isempty(infected_t)
				throw(ArgumentError("infected_t must not be empty"))
			end
			if region.I_surge_max !== nothing && region.surge_time + 1 > length(infected_t)
				throw(ArgumentError("surge_time $(region.surge_time) lies beyond the series (days 0 to $(length(infected_t) - 1))"))
			end

		#	Trajectory constraints: peak, early-infection size, and slope
			peak_ok = within_bound(maximum(infected_t), region.I_peak_max)
			surge_ok = region.I_surge_max === nothing ? nothing :
					   within_bound(infected_t[region.surge_time + 1], region.I_surge_max)
			max_slope = length(infected_t) > 1 ? maximum(diff(infected_t)) : 0.0
			slope_ok = within_bound(max_slope, region.S_max)

		#	Richards constraints: prevalence and final size
			fit_finite = all(isfinite, (K, r, ti, a))
			if region.P_max_thresh === nothing
				prev_ok = nothing
			else
				p_max = fit_finite ? peak_prevalence(K, r, a, rec_scale, rec_shape) : NaN
				prev_ok = within_bound(p_max, region.P_max_thresh)
			end
			final_ok = within_bound(fit_finite ? C0 + K : NaN, region.K_max)

		#	Conjunction of the constraints; a failed fit is never in target
			checks = (peak_ok, surge_ok, slope_ok, prev_ok, final_ok)
			all_ok = fit_finite && all(c -> c === nothing || c, checks)

		#	Assembling Result
			return (peak=peak_ok, surge=surge_ok, slope=slope_ok, prevalence=prev_ok,
					final_size=final_ok, in_target=all_ok)
	end
	@doc raw"""
	**Description**
	Evaluate each constraint of a `TargetRegion` for one simulated run and report which ones hold. `in_target` returns only the overall result.

	**Usage**
	```julia
	evaluate_target(region, infected_t, K, r, ti, a; C0=0.0, rec_scale=nothing, rec_shape=nothing)
	```

	**Arguments**
	- `region::TargetRegion`: The acceptable outcome region.
	- `infected_t::Vector{Float64}`: Proportion currently infected on each day, starting at day 0, one value per day.
	- `K, r, ti, a::Float64`: Fitted Richards parameters, with `K` measured beyond `C0`.
	- `C0::Float64`: Known initial offset, seeds/N (default `0.0`).
	- `rec_scale, rec_shape`: Weibull recovery parameters, required when `P_max_thresh` is set.

	**Details**
	The five constraints:
	1. **Peak infection size**: ``\max_t I(t) \le`` `I_peak_max`.
	2. **Early-infection size**: ``I(t_{\text{surge}}) \le`` `I_surge_max`, read at index `surge_time + 1`.
	3. **Slope**: ``\max_t \Delta I(t) \le`` `S_max`.
	4. **Prevalence**: `peak_prevalence(K, r, a, rec_scale, rec_shape)` ``\le`` `P_max_thresh`.
	5. **Final size**: ``C_0 + K \le`` `K_max`, the total proportion ever infected.

	Constraints 1 to 3 come from the trajectory; 4 and 5 come from the Richards fit. Constraint 1 is the observed peak; constraint 4 is the model-based peak from the fitted curve and the recovery distribution. If any Richards parameter is non-finite (a failed or extinct fit), `in_target` is `false`.

	`ArgumentError` is thrown for configuration errors: `P_max_thresh` set without the recovery parameters, an empty series, or a `surge_time` beyond the series.

	**Value**
	`NamedTuple` `(peak, surge, slope, prevalence, final_size, in_target)`. Each constraint is `true`, `false`, or `nothing` when unconstrained.

	**See Also**
	`in_target`, `TargetRegion`, `peak_prevalence`
	""" evaluate_target

#	Target Region Evaluation
	function in_target(region::TargetRegion,
					   infected_t::Vector{Float64},
					   K::Float64, r::Float64, ti::Float64, a::Float64;
					   C0::Float64=0.0,
					   rec_scale::Union{Float64, Nothing}=nothing,
					   rec_shape::Union{Float64, Nothing}=nothing)::Bool
		"""
		Args:
			Same as evaluate_target.
		Returns:
			Bool: true if every set constraint holds and the fit is finite
		Notes:
			Signature follows roadmap §4.3, with C0 added for Option A.
		"""

		#	Delegate to the detailed evaluation
			return evaluate_target(region, infected_t, K, r, ti, a;
								   C0=C0, rec_scale=rec_scale, rec_shape=rec_shape).in_target
	end
	@doc raw"""
	**Description**
	The binary indicator ``\mathbf{1}[\mathbf{x} \in \mathcal{A}]`` for one simulated run, stored as `in_target` in the output record.

	**Usage**
	`in_target(region, infected_t, K, r, ti, a; C0=0.0, rec_scale=nothing, rec_shape=nothing)`

	**Details**
	Returns `evaluate_target(...).in_target`. An all-`nothing` region returns `true` for any run with a finite fit.

	**Value**
	`Bool`.

	**See Also**
	`evaluate_target`, `TargetRegion`
	""" in_target

#   Exporting Objects
	export simulate_richards,
		   fit_richards,
		   richards_curve,
		   richards_peak,
		   TargetRegion,
		   discrete_weibull_mean,
		   peak_prevalence,
		   evaluate_target,
		   in_target

end # module

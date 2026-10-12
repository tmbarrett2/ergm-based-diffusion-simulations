#Richards Model Simulation Extension
#Jonathan H. Morgan, Ph.D.
#12 April 2026
#Revised 12 October 2026: known initial offset C(0) (Option A), convergence,
#boundary, and transition checks, fit status codes, peak timing and incidence,
#and an explicit random number generator for simulate_richards. Target region,
#analytic peak prevalence, and in_target implemented (UT-1.3). Non-starter
#runs (no transmission at all) are classified as non-viable; extinct runs
#follow an explicit policy in evaluate_target. The curve is parameterized by its
#inflection time t_i, the day of peak incidence, which stays within the horizon
#for early-peaking network epidemics. The fit also reports its root-mean-square
#residual, fit_rmse. The shape parameter a runs from -1 to 10, so the family
#spans the monomolecular (a = -1), Gompertz (a = 0), and logistic (a = 1) curves;
#a < 0 gives a curve with a definite start and an early peak. The share of the
#final size reached at the inflection is reported as frac_at_inflection.

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
			t_i::Float64: inflection time, the day of peak incidence
			a::Float64: asymmetry exponent (a ≥ -1; -1 monomolecular, 0 Gompertz, 1 logistic)
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
			C(t) = C0 + K [1 + a exp(-r (t - t_i))]^(-1/a)
			t_i is the inflection of C(t) and the day of peak incidence for
			every a > -1. Versions before 12 October 2026 took t_i as the
			location constant of the form K / [1 + exp(-r (t - t_i))]^(1/a),
			which is the inflection only when a = 1. For a ≠ 1, a curve
			generated with the same arguments is now shifted by ln(a)/r days.
			a may be 0 or negative down to -1 (12 October 2026); see
			richards_value for the forms these take.
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
			if !(a >= -1.0) || !isfinite(a)
				throw(DomainError(a, "a must be finite and at least -1"))
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
	- `t_i::Float64`: Inflection time, the day of peak incidence.
	- `a::Float64`: Asymmetry exponent, at least `-1`. `a = 1` is the logistic curve, `a = 0` the Gompertz curve, and `a = -1` the monomolecular curve.
	- `C0::Float64`: Known initial offset (default `0.0`). In simulator output this is the seed fraction, seeds/N.
	- `sigma::Union{Float64, Nothing}`: Standard deviation of additive Gaussian noise (default `nothing` for deterministic output).
	- `seed::Union{Int, Nothing}`: Seed for a local `Xoshiro` generator when noise is applied (default `nothing`).
	- `rng::Union{AbstractRNG, Nothing}`: Generator for the noise when `seed` is not set (default `nothing`, the global generator).

	**Details**
	The curve is:

	``C(t) = C_0 + K \left[1 + a\, e^{-r (t - t_i)}\right]^{-1/a}``

	Here ``t_i`` is the inflection of ``C(t)``, which is also the day of peak incidence, for every ``a > -1``. Before 12 October 2026, ``t_i`` was the location constant of the form ``K / [1 + e^{-r(t - t_i)}]^{1/a}``, which is the inflection only when ``a = 1``; for ``a \neq 1``, a curve generated with the same arguments is now shifted by ``\ln(a)/r`` days.

	Values ``-1 \le a \le 0`` are allowed; see `fit_richards` for what they mean.

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

#	Helper Function for richards_curve: the curve at one time point
	function richards_value(t::Real, K::Real, r::Real, ti::Real, a::Real)
		"""
		Args:
			t::Real: time point
			K, r, ti, a::Real: final size, growth rate, inflection time, and shape
		Returns:
			Float64 (or the promoted type): K [1 + a x]^(-1/a) with x = exp(-r (t - t_i))
		Notes:
			Three cases (12 October 2026):
			a ≠ 0: K exp(-log1p(a x)/a). log1p keeps small |a| accurate,
				so the curve is continuous through a = 0.
			a = 0: the Gompertz curve K exp(-x), the limit of the
				general form.
			a < 0 with a x ≤ -1: 0.0. A negative a gives a curve with a
				definite start at t_0 = t_i + ln(-a)/r (richards_start);
				before t_0 no one beyond the seeds is infected. The guard
				also keeps log1p inside its domain.
			Overflow of exp() gives x = Inf, and the result is 0.0 in
			every case, the correct limit.
		"""

		#	Gompertz case
			x = exp(-r * (t - ti))
			if a == 0
				return K * exp(-x)
			end

		#	Before the start of a curve with a definite start
			if a * x <= -1
				return zero(K * x)
			end

		#	General case
			return K * exp(-log1p(a * x) / a)
	end

#	Richards Curve
	function richards_curve(t, p)
		"""
		Args:
			t: time points (vector or scalar)
			p: parameter vector [K, r, t_i, a]
		Returns:
			Vector{Float64} (or Float64): K [1 + a exp(-r (t - t_i))]^(-1/a)
		Notes:
			Model function passed to LsqFit.curve_fit. t_i is the inflection
			time: setting the second derivative to zero gives
			a exp(-r (t - t_i)) = a, so the inflection, and the peak of
			incidence, falls at t = t_i for every a > -1. At a = -1 incidence
			is largest at the start, which is t_i.
			The shape a runs from -1 to 10 (12 October 2026). Evaluation is
			delegated to richards_value, which handles a = 0 (Gompertz) and
			the zero stretch before the start of a curve with a < 0.
		"""

		#	Unpack and evaluate point by point
			K, r, ti, a = p
			return richards_value.(t, K, r, ti, a)
	end
	@doc raw"""
	**Description**
	The Richards curve in its inflection-time form, the model fitted by `fit_richards`.

	**Usage**
	`richards_curve(t, [K, r, t_i, a])`

	**Details**
	``y(t) = K \left[1 + a\, e^{-r (t - t_i)}\right]^{-1/a}, \qquad -1 \le a \le 10.``

	Three classic growth curves are members of the family:

	| a | Curve | Formula | Incidence |
	|---|-------|---------|-----------|
	| 1 | Logistic | ``K / (1 + e^{-r(t - t_i)})`` | Symmetric; peak at half of ``K`` |
	| 0 | Gompertz | ``K \exp(-e^{-r(t - t_i)})`` | Peak at ``K/e``, about 37% of ``K`` |
	| -1 | Monomolecular | ``K \max(0, 1 - e^{-r(t - t_i)})`` | Largest at the start, ``t_i``, then declining |

	For ``a < 0`` the curve is exactly 0 until the start time ``t_0 = t_i + \ln(-a)/r`` (see `richards_start`) and rises like ``(t - t_0)^{-1/a}`` after it. For ``-1 < a`` the inflection, and so the peak of incidence, is at ``t_i``.

	**Value**
	`Float64`, or a vector of them when `t` is a vector.

	**See Also**
	`fit_richards`, `richards_start`, `inflection_fraction`
	""" richards_curve

#	Start Time of a Richards Curve
	function richards_start(r::Float64, t_i::Float64, a::Float64)::Float64
		"""
		Args:
			r::Float64: growth rate
			t_i::Float64: inflection time
			a::Float64: asymmetry exponent
		Returns:
			Float64: t_0 = t_i + ln(-a)/r for a < 0, or -Inf for a ≥ 0
		Notes:
			A curve with a ≥ 0 approaches 0 only as t → -∞, so it has no
			start; -Inf records that. For a < 0 the curve is 0 before t_0.
			Added 12 October 2026 with the negative range of a.
		"""

		#	Finite start only for a negative shape
			return a < 0.0 ? t_i + log(-a) / r : -Inf
	end
	@doc raw"""
	**Description**
	The day a Richards curve with ``a < 0`` begins to rise.

	**Usage**
	`richards_start(r, t_i, a)`

	**Details**
	For ``a < 0``, ``1 + a e^{-r(t - t_i)}`` reaches 0 at ``t_0 = t_i + \ln(-a)/r``, and the curve is 0 before that day. At ``a = -1``, ``t_0 = t_i``. For ``a \ge 0`` there is no finite start, and the function returns `-Inf`. For a simulated epidemic seeded on day 0, a fitted ``t_0`` near 0 is expected; a ``t_0`` well after 0 points to a smoldering phase before takeoff.

	**Value**
	`Float64`.

	**See Also**
	`richards_curve`, `fit_richards`
	""" richards_start

#	Share of the Final Size Reached at the Inflection
	function inflection_fraction(a::Float64)::Float64
		"""
		Args:
			a::Float64: asymmetry exponent (≥ -1)
		Returns:
			Float64: y(t_i)/K = (1 + a)^(-1/a)
		Notes:
			Computed as exp(-log1p(a)/a) for accuracy at small |a|.
			Equals 1/e at a = 0 (Gompertz), 0 at a = -1 (monomolecular),
			and 1/2 at a = 1 (logistic); it rises with a toward 1.
			Added 12 October 2026 as the shape summary for the run record.
		"""

		#	Gompertz and monomolecular limits, then the general case
			if a == 0.0
				return exp(-1.0)
			elseif a == -1.0
				return 0.0
			end
			return exp(-log1p(a) / a)
	end
	@doc raw"""
	**Description**
	The share of the final size reached by the inflection time, a direct reading of the curve's shape.

	**Usage**
	`inflection_fraction(a)`

	**Details**
	``\phi(a) = y(t_i)/K = (1 + a)^{-1/a}``, which increases with ``a``:

	| a | -1 | -0.5 | 0 | 1 | 10 |
	|---|----|------|---|---|----|
	| ``\phi`` | 0 | 0.25 | 0.37 | 0.50 | 0.79 |

	A small ``\phi`` means a front-loaded epidemic: incidence peaks early, before most infections occur, and is followed by a long tail. Unlike ``a``, ``\phi`` is bounded in ``[0, 1)`` and has a direct reading, so it is reported in the run record as `frac_at_inflection`.

	**Value**
	`Float64` in ``[0, 1)``.

	**See Also**
	`richards_curve`, `fit_richards`
	""" inflection_fraction

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

#	Richards Peak Incidence
	function peak_incidence(K::Float64, r::Float64, a::Float64)::Float64
		"""
		Args:
			K::Float64: final size (beyond the offset)
			r::Float64: growth rate
			a::Float64: asymmetry exponent
		Returns:
			Float64: the largest value of dC/dt, K r (1 + a)^(-(1 + a)/a)
		Notes:
			dC/dt peaks at the inflection time t_i, so its timing needs no
			separate function. The height does not depend on t_i.
			Replaces richards_peak, which returned the peak time from the
			former location constant; that time is now t_i itself.
			Computed as exp(-(1 + a) log1p(a)/a). The limits are taken
			exactly: K r/e at a = 0 (Gompertz) and K r at a = -1
			(monomolecular, where incidence is largest at the start).
		"""

		#	Gompertz and monomolecular limits
			if a == 0.0
				return K * r * exp(-1.0)
			elseif a == -1.0
				return K * r
			end

		#	Height of the incidence peak
			return K * r * exp(-(1.0 + a) * log1p(a) / a)
	end
	@doc raw"""
	**Description**
	Height of peak incidence for a Richards curve.

	**Usage**
	`peak_incidence(K, r, a)`

	**Details**
	For ``C(t) = K[1 + a e^{-r(t - t_i)}]^{-1/a}``, the incidence ``\dot C(t)`` is maximal at the inflection time ``t_i``, with

	``\dot C(t_i) = K r (1 + a)^{-(1 + a)/a}.``

	For ``a = 1`` this is ``Kr/4``; at ``a = 0`` (Gompertz) it is ``Kr/e``; at ``a = -1`` (monomolecular) it is ``Kr``, reached at the start of the curve. Under noise, ``r`` and ``a`` are individually weakly identified, but the peak incidence is recovered reliably, so it is reported with every fit.

	**Value**
	`Float64`.

	**See Also**
	`fit_richards`, `peak_prevalence`
	""" peak_incidence

#	Helper Function for fit_richards: assemble a fit result
	function fit_result(K::Float64, r::Float64, ti::Float64, a::Float64,
						fit_ok::Bool, extinct::Bool, status::Symbol; rmse::Float64=NaN)
		"""
		Args:
			K, r, ti, a::Float64: fitted (or NaN) Richards parameters
			fit_ok::Bool: true only for a usable fit
			extinct::Bool: true if the run was classified as extinct
			status::Symbol: reason code (see fit_richards)
			rmse::Float64: root-mean-square residual of a completed fit (default = NaN)
		Returns:
			NamedTuple: K, r_growth, t_i, a, frac_at_inflection, peak_incidence, fit_rmse,
						fit_ok, extinct, viable, status
		Notes:
			Replaces the nested result constructor that used to sit inside
			fit_richards. Peak incidence and the inflection fraction are
			computed whenever all four parameters are finite and a ≥ -1.
			t_i is the inflection time, so the former t_peak field would
			duplicate it and has been removed.
		"""

		#	Derived shape and peak quantities
			if all(isfinite, (K, r, ti, a)) && r > 0.0 && a >= -1.0
				frac = inflection_fraction(a)
				peak_inc = peak_incidence(K, r, a)
			else
				frac = NaN
				peak_inc = NaN
			end

		#	Assembling Result
			return (K=K, r_growth=r, t_i=ti, a=a, frac_at_inflection=frac, peak_incidence=peak_inc,
					fit_rmse=rmse, fit_ok=fit_ok, extinct=extinct, viable=(status != :nonstarter),
					status=status)
	end

#	Fit Richards Growth Model
	function fit_richards(t::Vector{Float64}, C::Vector{Float64};
						  C0::Float64=0.0,
						  extinction_threshold::Float64=0.02,
						  min_points::Int=5,
						  min_transition_points::Int=3,
						  bound_rtol::Float64=1e-6,
						  nonstarter_tol::Float64=1e-12,
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
			nonstarter_tol::Float64: C[end] - C0 at or below this → non-starter (default = 1e-12)
			max_iter::Int: iteration limit passed to LsqFit (default = 2000)
		Returns:
			NamedTuple with fields:
				K::Float64: fitted final size beyond C0
				r_growth::Float64: fitted per-capita growth rate
				t_i::Float64: fitted inflection time, the day of peak incidence
				a::Float64: fitted asymmetry exponent, in [-1, 10]
				frac_at_inflection::Float64: share of K reached at t_i, (1 + a)^(-1/a)
				peak_incidence::Float64: peak of dC/dt, K r (1 + a)^(-(1 + a)/a)
				fit_rmse::Float64: root-mean-square residual of the fit (NaN when no fit completed)
				fit_ok::Bool: true only for a usable fit
				extinct::Bool: true if 0 < C[end] - C0 < extinction_threshold
				viable::Bool: false for a non-starter (no transmission at all)
				status::Symbol: :ok, :nonstarter, :extinct, :too_few_points,
								:invalid_data, :no_convergence, :boundary,
								:too_steep, :implausible, or :error
		Notes:
			Option A (open decision 1, 12 October 2026): the known offset C0
			is subtracted before fitting, so K is the share infected beyond
			the seeds and extinction is judged on growth beyond the seeds.
			Conditioning assumption (12 October 2026): a run with no
			transmission at all is not an exercise of the model. It is
			classified as a non-starter (viable = false) and is excluded
			from fitting, training, and target evaluation.
			Inflection form (12 October 2026): the fitted parameters are
			(K, r, t_i, a), with t_i the inflection time. The former location
			constant, bounded below by 0, rejected 75 of 100 simulated village
			epidemics whose early peaks put that constant just below 0. The
			inflection time always lies within the horizon.
			Negative shapes (12 October 2026): with a > 0 the curve cannot
			inflect before reaching 1/e of K, and 64 of 85 established
			village runs pinned a at its former lower bound of 1e-6, with
			the largest residuals and inflection times later than the
			observed incidence peak. a now runs from -1 to 10. Gompertz
			(a = 0) is an interior point. Values below 0 describe a
			front-loaded epidemic: a definite start, an early peak, and a
			long tail. a = -1, the monomolecular curve, is a natural limit
			(no acceleration at all) and is not flagged; a fit near -1 is a
			result to examine, since it implies incidence driven by the
			remaining susceptibles alone.
			fit_rmse is computed over every day of the series, so a long flat
			tail after extinction lowers it; compare it across runs of the
			same horizon.
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

		#	Check for a non-starter: the seeds infected no one
			if y[end] <= nonstarter_tol
				return fit_result(NaN, NaN, NaN, NaN, false, false, :nonstarter)
			end

		#	Check for extinction: growth beyond the seeds is below threshold
			if y[end] < extinction_threshold
				return fit_result(NaN, NaN, NaN, NaN, false, true, :extinct)
			end

		#	Compute initial parameter guesses from data
		#	t_i0 uses the first crossing of half the final size (the inflection of a
		#	logistic curve, a0 = 1), which is less sensitive to noise than the
		#	largest single-day increase
			K0 = y[end]
			r0 = 0.3
			half_idx = findfirst(>=(K0 / 2.0), y)
			t_i0 = half_idx === nothing ? t[div(length(t), 2)] : t[half_idx]
			a0 = 1.0
			p0 = [K0, r0, t_i0, a0]

		#	Set parameter bounds on (K, r, t_i, a); K cannot exceed 1 - C0
		#	a = -1 (monomolecular) is the natural lower limit of the family
			T_max = t[end]
			lower = [0.0,      0.001, 0.0,   -1.0]
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

				#	Root-mean-square residual over the series
					rmse = sqrt(sum(abs2, fit.resid) / length(fit.resid))

				#	Check convergence: LsqFit returns rather than throws when it stops early
					if !fit.converged
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :no_convergence; rmse=rmse)
					end

				#	Check for implausible fit
					if K_fit < extinction_threshold
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :implausible; rmse=rmse)
					end

				#	Check bounds on r, t_i, and the upper bound of a
				#	K may reach 1 - C0 (everyone infected); a may reach -1 (monomolecular)
					if at_bound(r_fit, lower[2], upper[2]; rtol=bound_rtol) ||
					   at_bound(ti_fit, lower[3], upper[3]; rtol=bound_rtol) ||
					   (upper[4] - a_fit) <= bound_rtol * (upper[4] - lower[4])
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :boundary; rmse=rmse)
					end

				#	Check that the rise is resolved by the time grid
				#	A near-vertical step leaves no observations between 10% and 90% of K
					n_transition = count(v -> 0.1 * K_fit <= v <= 0.9 * K_fit, y)
					if n_transition < min_transition_points
						return fit_result(K_fit, r_fit, ti_fit, a_fit, false, false, :too_steep; rmse=rmse)
					end

				#	Assembling Result
					return fit_result(K_fit, r_fit, ti_fit, a_fit, true, false, :ok; rmse=rmse)

			catch e
				#	Numerical failure inside LsqFit
					return fit_result(NaN, NaN, NaN, NaN, false, false, :error)
			end
	end
	@doc raw"""
	**Description**
	Fit the Richards growth model to cumulative incidence data using nonlinear least squares, after subtracting the known initial offset. Returns fitted parameters, peak incidence, the residual error, and diagnostic flags for degenerate trajectories.

	**Usage**
	```julia
	fit_richards(t, C; C0=0.0, extinction_threshold=0.02, min_points=5,
				 min_transition_points=3, bound_rtol=1e-6, nonstarter_tol=1e-12,
				 max_iter=2000)
	```

	**Arguments**
	- `t::Vector{Float64}`: Time points of the observations.
	- `C::Vector{Float64}`: Cumulative proportion ever infected at each time point.
	- `C0::Float64`: Known initial offset (default `0.0`). For simulator output, pass the seed fraction, seeds/N.
	- `extinction_threshold::Float64`: If `C[end] - C0` falls below this value, the run is classified as extinct and fitting is skipped (default `0.02`).
	- `min_points::Int`: Minimum number of data points required for fitting (default `5`).
	- `min_transition_points::Int`: Minimum number of observations between 10% and 90% of the fitted `K` (default `3`).
	- `bound_rtol::Float64`: Tolerance for the boundary check, as a fraction of each parameter's bound range (default `1e-6`).
	- `nonstarter_tol::Float64`: If `C[end] - C0` is at or below this value, the seeds infected no one and the run is a non-starter (default `1e-12`).
	- `max_iter::Int`: Iteration limit passed to `LsqFit` (default `2000`).

	**Details**
	The model fitted to ``y(t) = C(t) - C_0`` is the Richards curve written in terms of its inflection time:

	``y(t) = K \left[1 + a\, e^{-r (t - t_i)}\right]^{-1/a}``

	The inflection of ``y``, and so the peak of incidence, falls at ``t_i`` for every ``a > -1``. For ``a > 0`` the curve is identical to ``K / [1 + e^{-r(t - \ell)}]^{1/a}`` with location constant ``\ell = t_i + \ln(a)/r``, which is the inflection only when ``a = 1``. For early-peaking epidemics ``\ell`` falls below 0, so it is not reported.

	The shape ``a`` runs from ``-1`` to ``10``, so the family includes the logistic (``a = 1``), Gompertz (``a = 0``), and monomolecular (``a = -1``) curves:

	| a | Share of ``K`` at the peak of incidence | Reading |
	|---|------------------------------------------|---------|
	| 1 | 0.50 | Symmetric rise and fall, as under homogeneous mixing |
	| 0 | 0.37 | Brief acceleration, then a long tail |
	| -0.5 | 0.25 | Front-loaded: a definite start, an early peak, a longer tail |
	| -1 | 0 | No acceleration: incidence is largest on the first day |

	A curve with ``a < 0`` is 0 until ``t_0 = t_i + \ln(-a)/r`` (`richards_start`), so it has a definite start, as a single introduced seed does. Before 12 October 2026 the lower bound on ``a`` was ``10^{-6}``; with ``a > 0`` the curve cannot inflect before reaching ``1/e`` of ``K``, and most simulated village epidemics pinned ``a`` at that bound, with larger residuals and a late ``t_i``. Front-loading changes timing and asymmetry far more than height: at a fixed final size and rise time, the peak incidence at ``a = -0.5`` is only about 18% above the logistic value.

	The share of ``K`` reached at the inflection, ``(1 + a)^{-1/a}``, is reported as `frac_at_inflection` (`inflection_fraction`). It is bounded in ``[0, 1)`` and reads directly, so it is the preferred shape summary.

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
	| t_i | 0.0 | T_max | Inflection time within the horizon |
	| a | -1.0 | 10.0 | Shape: 1 logistic, 0 Gompertz, -1 monomolecular |

	Degenerate trajectory handling (`status` in parentheses):
	- **Too few points** (`:too_few_points`): `length(t) < min_points`.
	- **Non-starter** (`:nonstarter`): `C[end] - C0 ≤ nonstarter_tol`, so the seeds infected no one; `viable=false`. Every analysis is conditioned on at least one transmission, so a non-starter is not an exercise of the model and is excluded throughout.
	- **Invalid data** (`:invalid_data`): mismatched lengths, non-finite values, or `C0` outside `[0, 1)`.
	- **Extinction** (`:extinct`): `C[end] - C0 < extinction_threshold`; `extinct=true`.
	- **Non-convergence** (`:no_convergence`): `LsqFit` stopped without meeting its convergence tolerances. `LsqFit` reports this through `converged = false` rather than by throwing.
	- **Implausible fit** (`:implausible`): fitted `K < extinction_threshold`.
	- **Boundary fit** (`:boundary`): `r` or `t_i` on either bound, or `a` on its upper bound. `K` may reach `1 - C0`, which means everyone beyond the seeds was infected, and `a` may reach `-1`, the monomolecular curve. Both are natural limits of what they measure, so neither is flagged.
	- **Too steep** (`:too_steep`): fewer than `min_transition_points` observations between 10% and 90% of the fitted `K`, so the rise is not resolved by the time grid.
	- **Numerical error** (`:error`): `LsqFit` threw, or returned non-finite parameters.

	No exception propagates out of `fit_richards`. When `fit_ok=false`, parameters are `NaN` for failures that occur before or during fitting, and the fitted values are returned for diagnostic inspection when a fit completed but failed a check.

	Under noise, `r` and `a` trade off against one another and are individually weakly identified. `K`, `t_i`, and `peak_incidence` are recovered reliably.

	`fit_rmse` is the root-mean-square difference between the fitted curve and the data over every day of the series. It is reported for every fit that completed, including fits that fail a later check, and is `NaN` otherwise. Because the days after extinction are flat and well fitted, it should be compared across runs of the same horizon.

	**Value**
	A `NamedTuple` with fields:
	- `K::Float64`: Fitted final size beyond `C0`.
	- `r_growth::Float64`: Fitted per-capita growth rate.
	- `t_i::Float64`: Fitted inflection time, the day of peak incidence.
	- `a::Float64`: Fitted asymmetry exponent, in ``[-1, 10]``.
	- `frac_at_inflection::Float64`: Share of `K` reached at `t_i`, ``(1 + a)^{-1/a}``.
	- `peak_incidence::Float64`: Peak of ``\dot C``.
	- `fit_rmse::Float64`: Root-mean-square residual of the fit.
	- `fit_ok::Bool`: `true` if fitting succeeded and passed every check.
	- `extinct::Bool`: `true` if the run was classified as extinct (some transmission, but growth beyond the seeds below `extinction_threshold`).
	- `viable::Bool`: `false` for a non-starter.
	- `status::Symbol`: Reason code, `:ok` for a usable fit.

	**Examples**
	```julia
	#	Generate a known curve with three seeds on 150 nodes
		t = collect(0.0:1.0:200.0)
		result = simulate_richards(t, 0.40, 0.15, 70.0, 1.0; C0=3/150)

	#	Recover parameters, passing the known offset
		fit = fit_richards(t, result.C; C0=3/150)
		println("K: ", fit.K, "  r: ", fit.r_growth, "  t_i: ", fit.t_i, "  a: ", fit.a)
		println("fit_rmse: ", fit.fit_rmse, "  status: ", fit.status)
	```

	**References**
	Richards, F.J. (1959). A flexible growth function for empirical use.
	*Journal of Experimental Botany*, 10(2), 290–301.

	**See Also**
	`simulate_richards`, `richards_curve`, `richards_start`, `inflection_fraction`, `peak_incidence`, `LsqFit.curve_fit`
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
			Richards curve (peak_incidence; t_i does not affect its height).
			E[D] is the mean of the discretized recovery duration
			(discrete_weibull_mean). P_max approximates peak prevalence by
			treating incidence as constant near the peak for one mean
			infectious period; it overstates the true peak when the
			incidence peak is narrow relative to E[D].
		"""

		#	Peak incidence times mean infectious duration
			I_max = peak_incidence(K, r, a)
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
	`peak_incidence`, `discrete_weibull_mean`, `in_target`
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
							 rec_shape::Union{Float64, Nothing}=nothing,
							 status::Symbol=:ok,
							 extinct_policy::Symbol=:exclude,
							 final_cum::Union{Float64, Nothing}=nothing)
		"""
		Args:
			region::TargetRegion: acceptable outcome region
			infected_t::Vector{Float64}: proportion currently infected on days 0, 1, ..., T
										 (one value per day; summarize.jl expands the simulator log)
			K, r, ti, a::Float64: fitted Richards parameters (K beyond C0; ti the inflection time; NaN if not fitted)
			C0::Float64: known initial offset, seeds/N (default = 0.0)
			rec_scale, rec_shape: Weibull recovery parameters, required when P_max_thresh is set
			status::Symbol: status from fit_richards (default = :ok)
			extinct_policy::Symbol: :exclude, :success, or :failure for extinct runs (default = :exclude)
			final_cum::Union{Float64, Nothing}: C(end), the final cumulative proportion;
												needed to judge final size for an extinct run
		Returns:
			NamedTuple: (peak, surge, slope, prevalence, final_size, in_target, class)
			Each constraint is true, false, or nothing (unconstrained or not evaluated).
			in_target is true, false, or nothing (excluded). class is :established,
			:extinct, :nonstarter, or :fit_failed.
		Notes:
			Every analysis is conditioned on at least one transmission
			(12 October 2026). A non-starter is not an exercise of the model,
			so in_target is nothing under every policy. Extinct runs follow
			extinct_policy: :exclude (default) conditions on establishment,
			:success counts a contained run as acceptable when it meets the
			trajectory and final-size bounds, and :failure reproduces the
			original roadmap rule. A failed fit of an established run is
			excluded (nothing), because the failure reflects the curve fit
			rather than the behavior being modeled; fit failure rates are
			reported separately.
			Throws ArgumentError for configuration errors: an unknown policy,
			P_max_thresh without recovery parameters, an empty series, a
			surge_time beyond the series, or an extinct run judged under
			:success with K_max set but no final_cum.
		"""

		#	Validation of the configuration
			if !(extinct_policy in (:exclude, :success, :failure))
				throw(ArgumentError("extinct_policy must be :exclude, :success, or :failure"))
			end
			if region.P_max_thresh !== nothing && (rec_scale === nothing || rec_shape === nothing)
				throw(ArgumentError("P_max_thresh is set, so rec_scale and rec_shape must be provided"))
			end
			if isempty(infected_t)
				throw(ArgumentError("infected_t must not be empty"))
			end
			if region.I_surge_max !== nothing && region.surge_time + 1 > length(infected_t)
				throw(ArgumentError("surge_time $(region.surge_time) lies beyond the series (days 0 to $(length(infected_t) - 1))"))
			end

		#	Classify the run
			fit_finite = all(isfinite, (K, r, ti, a))
			if status == :nonstarter
				class = :nonstarter
			elseif status == :extinct
				class = :extinct
			elseif status == :ok && fit_finite
				class = :established
			else
				class = :fit_failed
			end

		#	Early return: non-starters and failed fits are not evaluated
			if class == :nonstarter || class == :fit_failed
				return (peak=nothing, surge=nothing, slope=nothing, prevalence=nothing,
						final_size=nothing, in_target=nothing, class=class)
			end

		#	Trajectory constraints: peak, early-infection size, and slope
			peak_ok = within_bound(maximum(infected_t), region.I_peak_max)
			surge_ok = region.I_surge_max === nothing ? nothing :
					   within_bound(infected_t[region.surge_time + 1], region.I_surge_max)
			max_slope = length(infected_t) > 1 ? maximum(diff(infected_t)) : 0.0
			slope_ok = within_bound(max_slope, region.S_max)

		#	Fit-based constraints for an established run
			if class == :established
				#	Prevalence from the fitted curve and the recovery distribution
					prev_ok = region.P_max_thresh === nothing ? nothing :
							  within_bound(peak_prevalence(K, r, a, rec_scale, rec_shape), region.P_max_thresh)

				#	Final size: total proportion ever infected
					final_ok = within_bound(C0 + K, region.K_max)

				#	Conjunction of the set constraints
					checks = (peak_ok, surge_ok, slope_ok, prev_ok, final_ok)
					all_ok = all(c -> c === nothing || c, checks)
					return (peak=peak_ok, surge=surge_ok, slope=slope_ok, prevalence=prev_ok,
							final_size=final_ok, in_target=all_ok, class=class)
			end

		#	Extinct run: no fitted curve, so use the observed series
		#	Prevalence falls back to the observed peak; final size uses C(end)
			if extinct_policy == :success && region.K_max !== nothing && final_cum === nothing
				throw(ArgumentError("final_cum is required to judge final size for an extinct run under :success"))
			end
			prev_ok = region.P_max_thresh === nothing ? nothing :
					  within_bound(maximum(infected_t), region.P_max_thresh)
			final_ok = (region.K_max === nothing || final_cum === nothing) ? nothing :
					   within_bound(final_cum, region.K_max)

		#	Apply the extinct-run policy
			if extinct_policy == :exclude
				in_tgt = nothing
			elseif extinct_policy == :failure
				in_tgt = false
			else
				checks = (peak_ok, surge_ok, slope_ok, prev_ok, final_ok)
				in_tgt = all(c -> c === nothing || c, checks)
			end

		#	Assembling Result
			return (peak=peak_ok, surge=surge_ok, slope=slope_ok, prevalence=prev_ok,
					final_size=final_ok, in_target=in_tgt, class=class)
	end
	@doc raw"""
	**Description**
	Evaluate each constraint of a `TargetRegion` for one simulated run, classify the run, and report whether it lies in the acceptable region. `in_target` returns only the overall result.

	**Usage**
	```julia
	evaluate_target(region, infected_t, K, r, ti, a; C0=0.0, rec_scale=nothing, rec_shape=nothing,
					status=:ok, extinct_policy=:exclude, final_cum=nothing)
	```

	**Arguments**
	- `region::TargetRegion`: The acceptable outcome region.
	- `infected_t::Vector{Float64}`: Proportion currently infected on each day, starting at day 0, one value per day.
	- `K, r, ti, a::Float64`: Fitted Richards parameters, with `K` measured beyond `C0` and `ti` the inflection time (`NaN` when no fit exists).
	- `C0::Float64`: Known initial offset, seeds/N (default `0.0`).
	- `rec_scale, rec_shape`: Weibull recovery parameters, required when `P_max_thresh` is set.
	- `status::Symbol`: The `status` returned by `fit_richards` (default `:ok`).
	- `extinct_policy::Symbol`: `:exclude` (default), `:success`, or `:failure` for extinct runs.
	- `final_cum`: Final cumulative proportion C(end); required to judge final size for an extinct run under `:success`.

	**Details**
	Every analysis is conditioned on at least one transmission. Runs fall into four classes:

	| Class | Source | `in_target` |
	|-------|--------|-------------|
	| `:nonstarter` | `status == :nonstarter` (the seeds infected no one) | `nothing` under every policy |
	| `:extinct` | `status == :extinct` (some transmission, growth below the extinction threshold) | set by `extinct_policy` |
	| `:established` | `status == :ok` with finite parameters | evaluated on all set constraints |
	| `:fit_failed` | any other status | `nothing` (excluded) |

	For an established run, the five constraints are:
	1. **Peak infection size**: ``\max_t I(t) \le`` `I_peak_max`.
	2. **Early-infection size**: ``I(t_{\text{surge}}) \le`` `I_surge_max`, read at index `surge_time + 1`.
	3. **Slope**: ``\max_t \Delta I(t) \le`` `S_max`.
	4. **Prevalence**: `peak_prevalence(K, r, a, rec_scale, rec_shape)` ``\le`` `P_max_thresh`.
	5. **Final size**: ``C_0 + K \le`` `K_max`, the total proportion ever infected.

	An extinct run has no fitted curve, so constraint 4 uses the observed peak and constraint 5 uses `final_cum`. Under `:exclude`, which conditions on establishment, `in_target` is `nothing`; under `:failure` it is `false` (the original roadmap rule); under `:success` it is evaluated on those observed quantities.

	`ArgumentError` is thrown for configuration errors: an unknown policy, `P_max_thresh` set without the recovery parameters, an empty series, a `surge_time` beyond the series, or an extinct run judged under `:success` with `K_max` set but no `final_cum`.

	**Value**
	`NamedTuple` `(peak, surge, slope, prevalence, final_size, in_target, class)`. Each constraint is `true`, `false`, or `nothing`; `in_target` is `true`, `false`, or `nothing` when the run is excluded.

	**See Also**
	`in_target`, `TargetRegion`, `peak_prevalence`, `fit_richards`
	""" evaluate_target

#	Target Region Evaluation
	function in_target(region::TargetRegion,
					   infected_t::Vector{Float64},
					   K::Float64, r::Float64, ti::Float64, a::Float64;
					   C0::Float64=0.0,
					   rec_scale::Union{Float64, Nothing}=nothing,
					   rec_shape::Union{Float64, Nothing}=nothing,
					   status::Symbol=:ok,
					   extinct_policy::Symbol=:exclude,
					   final_cum::Union{Float64, Nothing}=nothing)::Union{Bool, Nothing}
		"""
		Args:
			Same as evaluate_target.
		Returns:
			Union{Bool, Nothing}: true or false for an evaluated run; nothing
			for a non-starter, a failed fit, or an extinct run under :exclude
		Notes:
			Signature follows roadmap §4.3, with C0 added for Option A and
			status, extinct_policy, and final_cum added for run classification.
			In the JSONL record, nothing is written as null.
		"""

		#	Delegate to the detailed evaluation
			return evaluate_target(region, infected_t, K, r, ti, a;
								   C0=C0, rec_scale=rec_scale, rec_shape=rec_shape,
								   status=status, extinct_policy=extinct_policy,
								   final_cum=final_cum).in_target
	end
	@doc raw"""
	**Description**
	The indicator ``\mathbf{1}[\mathbf{x} \in \mathcal{A}]`` for one simulated run, stored as `in_target` in the output record, or `nothing` (JSON `null`) when the run is excluded.

	**Usage**
	`in_target(region, infected_t, K, r, ti, a; C0=0.0, rec_scale=nothing, rec_shape=nothing, status=:ok, extinct_policy=:exclude, final_cum=nothing)`

	**Details**
	Returns `evaluate_target(...).in_target`. An all-`nothing` region returns `true` for any established run.

	**Value**
	`Union{Bool, Nothing}`.

	**See Also**
	`evaluate_target`, `TargetRegion`
	""" in_target

#   Exporting Objects
	export simulate_richards,
		   fit_richards,
		   richards_curve,
		   richards_start,
		   inflection_fraction,
		   peak_incidence,
		   TargetRegion,
		   discrete_weibull_mean,
		   peak_prevalence,
		   evaluate_target,
		   in_target

end # module

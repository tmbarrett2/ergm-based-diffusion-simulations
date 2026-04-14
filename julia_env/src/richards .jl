#Richards Model Simulation Extension
#Jonathan H. Morgan, Ph.D.
#12 April 2026

module richards_model

using LsqFit
using Distributions

# ============================================================
#	SIMULATION
# ============================================================

#	Simulate Richards Growth Curve
	function simulate_richards(t::Vector{Float64}, K::Float64, r::Float64, t_i::Float64, a::Float64;
	                           sigma::Union{Float64, Nothing}=nothing, seed::Union{Int, Nothing}=nothing)
		"""
		Args:
			t::Vector{Float64}: time points at which to evaluate the curve
			K::Float64: final epidemic size (asymptotic proportion, 0 < K ≤ 1)
			r::Float64: per-capita growth rate (r > 0)
			t_i::Float64: inflection time
			a::Float64: asymmetry exponent (a > 0; a = 1 recovers logistic)
			sigma::Union{Float64, Nothing}: standard deviation of additive Gaussian noise (default = nothing, deterministic)
			seed::Union{Int, Nothing}: random seed for reproducibility when sigma is set (default = nothing)
		Returns:
			NamedTuple with fields:
				t::Vector{Float64}: input time vector
				C::Vector{Float64}: cumulative incidence values (noisy if sigma is set)
				C_true::Vector{Float64}: deterministic Richards curve (always clean)
				params::NamedTuple: (K, r, t_i, a) used to generate the curve
		Notes:
			C(t) = K / [1 + exp(-r * (t - t_i))]^(1/a)
			When sigma is provided, Gaussian noise N(0, sigma^2) is added and
			values are clamped to [0.0, 1.0] since C represents a proportion.
		"""

		#	Validation
			if K <= 0.0 || K > 1.0
				throw(DomainError(K, "K must be in (0, 1]"))
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

		#	Compute deterministic Richards curve
			C_true = K ./ (1.0 .+ exp.(-r .* (t .- t_i))) .^ (1.0 / a)

		#	Add optional Gaussian noise
			if sigma !== nothing && sigma > 0.0
				#	Set random seed if provided
					if seed !== nothing
						Random.seed!(seed)
					end

				#	Generate noisy curve clamped to [0, 1]
					noise = sigma .* randn(length(t))
					C = clamp.(C_true .+ noise, 0.0, 1.0)
			else
				#	Deterministic case
					C = copy(C_true)
			end

		#	Assembling Result
			return (t=t, C=C, C_true=C_true, params=(K=K, r=r, t_i=t_i, a=a))
	end
	@doc raw"""
	**Description**
	Generate a Richards growth curve with known parameters, optionally perturbed by additive Gaussian noise. Used to produce synthetic cumulative incidence trajectories for testing `fit_richards` parameter recovery.

	**Usage**
	```julia
	simulate_richards(t, K, r, t_i, a; sigma=nothing, seed=nothing)
	```

	**Arguments**
	- `t::Vector{Float64}`: Time points at which to evaluate the curve.
	- `K::Float64`: Final epidemic size as a proportion (0 < K ≤ 1).
	- `r::Float64`: Per-capita growth rate (must be positive).
	- `t_i::Float64`: Inflection time of the sigmoid.
	- `a::Float64`: Asymmetry exponent (must be positive). `a = 1` recovers the standard logistic curve.
	- `sigma::Union{Float64, Nothing}`: Standard deviation of additive Gaussian noise (default `nothing` for deterministic output).
	- `seed::Union{Int, Nothing}`: Random seed for reproducibility when noise is applied (default `nothing`).

	**Details**
	The Richards growth model is defined as:

	``C(t) = \frac{K}{[1 + \exp(-r \cdot (t - t_i))]^{1/a}}``

	When `sigma` is `nothing` or `0.0`, the output is the exact deterministic curve. When `sigma > 0`, independent Gaussian noise ``N(0, \sigma^2)`` is added to each point and the result is clamped to `[0.0, 1.0]` since `C` represents a proportion.

	The deterministic curve is always returned in `C_true` regardless of the noise setting, enabling direct comparison between noisy observations and ground truth.

	**Value**
	A `NamedTuple` with fields:
	- `t::Vector{Float64}`: The input time vector.
	- `C::Vector{Float64}`: Cumulative incidence values (noisy if `sigma` is set, otherwise identical to `C_true`).
	- `C_true::Vector{Float64}`: The clean deterministic Richards curve.
	- `params::NamedTuple`: The generating parameters `(K, r, t_i, a)`.

	**Examples**
	```julia
	using Random

	#	Define time grid
		t = collect(0.0:1.0:200.0)

	#	Deterministic curve
		result = simulate_richards(t, 0.15, 0.10, 70.0, 1.3)
		println("Final size: ", result.C[end])

	#	Noisy curve with fixed seed
		noisy = simulate_richards(t, 0.15, 0.10, 70.0, 1.3; sigma=0.01, seed=42)
		println("Max deviation: ", maximum(abs.(noisy.C .- noisy.C_true)))
	```

	**References**
	Richards, F.J. (1959). A flexible growth function for empirical use.
	*Journal of Experimental Botany*, 10(2), 290–301.
	""" simulate_richards

# ============================================================
#	FITTING
# ============================================================

#	Fit Richards Growth Model
	function fit_richards(t::Vector{Float64}, C::Vector{Float64};
	                      extinction_threshold::Float64=0.02,
	                      min_points::Int=5)
		"""
		Args:
			t::Vector{Float64}: time points
			C::Vector{Float64}: cumulative proportion ever infected
			extinction_threshold::Float64: C[end] below this → extinct (default = 0.02)
			min_points::Int: minimum number of data points required (default = 5)
		Returns:
			NamedTuple with fields:
				K::Float64: fitted final epidemic size
				r_growth::Float64: fitted per-capita growth rate
				t_i::Float64: fitted inflection time
				a::Float64: fitted asymmetry exponent
				fit_ok::Bool: false if fitting failed for any reason
				extinct::Bool: true if C[end] < extinction_threshold
		Notes:
			Uses LsqFit.curve_fit with bounded parameters.
			Degenerate trajectories (extinction, too few points,
			non-convergence, implausible K) return fit_ok=false
			with NaN parameter values.
		"""

		#	Define failed result constructor
			failed = (K=NaN, r_growth=NaN, t_i=NaN, a=NaN, fit_ok=false, extinct=false)

		#	Check minimum data points
			if length(t) < min_points
				return failed
			end

		#	Check for extinction
			if C[end] < extinction_threshold
				return (K=NaN, r_growth=NaN, t_i=NaN, a=NaN, fit_ok=false, extinct=true)
			end

		#	Define Richards model for LsqFit
			function richards_model(t, p)
				K, r, ti, a = p
				return K ./ (1.0 .+ exp.(-r .* (t .- ti))) .^ (1.0 / a)
			end

		#	Compute initial parameter guesses from data
			K0 = C[end]
			r0 = 0.3
			dC = diff(C)
			t_i0 = length(dC) > 0 ? t[argmax(dC)] : t[div(length(t), 2)]
			a0 = 1.0
			p0 = [K0, r0, t_i0, a0]

		#	Set parameter bounds
			T_max = t[end]
			lower = [0.0,   0.001, 0.0,  0.1]
			upper = [1.0,   5.0,   T_max, 10.0]

		#	Clamp initial guesses to bounds
			p0 = clamp.(p0, lower, upper)

		#	Attempt curve fitting
			try
				#	Run nonlinear least squares
					fit = curve_fit(richards_model, t, C, p0; lower=lower, upper=upper)

				#	Extract fitted parameters
					K_fit = fit.param[1]
					r_fit = fit.param[2]
					ti_fit = fit.param[3]
					a_fit = fit.param[4]

				#	Check for implausible fit
					if K_fit < extinction_threshold
						return (K=K_fit, r_growth=r_fit, t_i=ti_fit, a=a_fit, fit_ok=false, extinct=false)
					end

				#	Assembling Result
					return (K=K_fit, r_growth=r_fit, t_i=ti_fit, a=a_fit, fit_ok=true, extinct=false)

			catch e
				#	Non-convergence or numerical failure
					return failed
			end
	end
	@doc raw"""
	**Description**
	Fit the Richards growth model to cumulative incidence data using nonlinear least squares. Returns fitted parameters and diagnostic flags for degenerate trajectories (extinction, insufficient data, non-convergence, implausible fits).

	**Usage**
	```julia
	fit_richards(t, C; extinction_threshold=0.02, min_points=5)
	```

	**Arguments**
	- `t::Vector{Float64}`: Time points of the observations.
	- `C::Vector{Float64}`: Cumulative proportion ever infected at each time point.
	- `extinction_threshold::Float64`: If `C[end]` falls below this value, the epidemic is classified as extinct and fitting is skipped (default `0.02`).
	- `min_points::Int`: Minimum number of data points required for fitting (default `5`).

	**Details**
	The Richards growth model is:

	``C(t) = \frac{K}{[1 + \exp(-r \cdot (t - t_i))]^{1/a}}``

	Initial parameter guesses are derived from the data:
	- ``K_0 = C[\text{end}]``
	- ``r_0 = 0.3``
	- ``t_{i,0} = t[\text{argmax}(\Delta C)]`` (time of steepest increase)
	- ``a_0 = 1.0``

	Parameter bounds enforced during fitting:

	| Parameter | Lower | Upper | Notes |
	|-----------|-------|-------|-------|
	| K | 0.0 | 1.0 | Proportion; hard physical bound |
	| r | 0.001 | 5.0 | Growth rate; upper bound prevents divergence |
	| t_i | 0.0 | T_max | Inflection within observation horizon |
	| a | 0.1 | 10.0 | Shape parameter; a=1 is logistic |

	Degenerate trajectory handling:
	- **Extinction**: `C[end] < extinction_threshold` → `extinct=true`, `fit_ok=false`
	- **Too few points**: `length(t) < min_points` → `fit_ok=false`
	- **Non-convergence**: `LsqFit` throws an exception → `fit_ok=false`
	- **Implausible fit**: fitted `K < extinction_threshold` → `fit_ok=false`

	When `fit_ok=false`, parameter values are `NaN` (except for the implausible-K case, which returns the fitted values for diagnostic inspection).

	**Value**
	A `NamedTuple` with fields:
	- `K::Float64`: Fitted final epidemic size.
	- `r_growth::Float64`: Fitted per-capita growth rate.
	- `t_i::Float64`: Fitted inflection time.
	- `a::Float64`: Fitted asymmetry exponent.
	- `fit_ok::Bool`: `true` if fitting succeeded and produced plausible parameters.
	- `extinct::Bool`: `true` if the trajectory was classified as extinct.

	**Examples**
	```julia
	using LsqFit, Random

	#	Generate a known curve
		t = collect(0.0:1.0:200.0)
		result = simulate_richards(t, 0.40, 0.15, 70.0, 1.0)

	#	Recover parameters
		fit = fit_richards(t, result.C)
		println("K: ", fit.K, "  r: ", fit.r_growth, "  t_i: ", fit.t_i, "  a: ", fit.a)
		println("fit_ok: ", fit.fit_ok, "  extinct: ", fit.extinct)
	```

	**References**
	Richards, F.J. (1959). A flexible growth function for empirical use.
	*Journal of Experimental Botany*, 10(2), 290–301.

	**See Also**
	`simulate_richards`, `LsqFit.curve_fit`
	""" fit_richards

# ============================================================
#	TARGET REGION
# ============================================================

#	Target Region Struct
	# TODO: Implement five epidemiological constraints:
	#   1. Peak infection size:    max(infected_t) ≤ I_peak_max
	#   2. Early-infection size:   infected[surge_time] ≤ I_surge_max
	#   3. Slope:                  max(Δ infected_t) ≤ S_max
	#   4. Prevalence:             P_max ≤ P_max_thresh
	#   5. Final size:             K ≤ K_max
	# All bounds optional; omitted = unconstrained.
	# surge_time defaults to 21 days.
	# struct TargetRegion
	#     I_peak_max   :: Union{Float64, Nothing}
	#     I_surge_max  :: Union{Float64, Nothing}
	#     surge_time   :: Int
	#     S_max        :: Union{Float64, Nothing}
	#     P_max_thresh :: Union{Float64, Nothing}
	#     K_max        :: Union{Float64, Nothing}
	# end

#	Analytic Peak Prevalence
	# TODO: Implement analytic peak prevalence from Richards
	# parameters and Weibull recovery distribution.
	# I_max = (K * r / a) * a^(a/(a+1)) / (1+a)^((a+1)/a)
	# P_max = I_max * E[D]
	# where E[D] = rec_scale * Γ(1 + 1/rec_shape)
	# For standard logistic (a=1): I_max = K*r/4
	# function peak_prevalence(K, r, a, rec_scale, rec_shape)::Float64

#	Target Region Evaluation
	# TODO: Implement in_target evaluation against five constraints.
	# Takes trajectory (infected_t), Richards params (K, r, ti, a),
	# and optional Weibull recovery params (rec_scale, rec_shape).
	# Throws ArgumentError if P_max_thresh is set but Weibull
	# params are not provided.
	# function in_target(region::TargetRegion,
	#                    infected_t::Vector{Float64},
	#                    K::Float64, r::Float64, ti::Float64, a::Float64;
	#                    rec_scale::Union{Float64,Nothing}=nothing,
	#                    rec_shape::Union{Float64,Nothing}=nothing)::Bool

end # module
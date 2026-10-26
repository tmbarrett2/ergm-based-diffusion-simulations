#Run Summary Pipeline: Daily Expansion, Richards Fit, and Target Evaluation
#Jonathan H. Morgan, Ph.D.
#12 October 2026
#Revised 12 October 2026: t_i is the inflection time (the day of peak incidence),
#so the separate t_peak field is removed; fit_rmse is added to the run summary.
#The shape a now runs from -1 to 10, and frac_at_inflection, the share of the
#final size reached at t_i, is added as the shape summary.
#Revised 26 October 2026: t_10, t_half, and t_90, the days on which the run
#reaches 10%, 50%, and 90% of its growth beyond the seeds, are taken from the
#data, because the fitted t_i marks the takeoff rather than the peak when a
#reaches its lower limit of -1.

module run_summary

using DataFrames
using ..richards_model

# ============================================================
#	DAILY EXPANSION
# ============================================================

#	Expand a Simulator Log to One Row per Day
	function expand_daily_log(ilog::DataFrame, maxtime::Int)
		"""
		Args:
			ilog::DataFrame: infection_log returned by sirdif
			maxtime::Int: simulation horizon in days
		Returns:
			DataFrame: one row for each day 0, 1, ..., maxtime with the same columns
		Notes:
			sirdif follows the SAS model: when the last infection clears it
			writes one row at maxtime and stops, so days between the last
			simulated day and maxtime are missing. Each missing day takes the
			values of the latest row at or before it. Nothing changes after
			extinction (no one is infected, so no one can be infected), so
			carrying the last state forward is exact for the columns the
			summaries use. The row at maxtime is kept as written.
			Throws ArgumentError if the log is empty, does not start at day 0,
			has non-integer or decreasing times, or runs past maxtime.
		"""

		#	Validation
			if nrow(ilog) == 0
				throw(ArgumentError("infection log is empty"))
			end
			if maxtime < 0
				throw(ArgumentError("maxtime must be ≥ 0"))
			end
			times = ilog.time
			if times[1] != 0.0
				throw(ArgumentError("infection log must start at day 0"))
			end
			if any(t -> t != round(t), times)
				throw(ArgumentError("infection log times must be whole days"))
			end
			if any(diff(times) .<= 0.0)
				throw(ArgumentError("infection log times must be strictly increasing"))
			end
			if times[end] > maxtime
				throw(ArgumentError("infection log runs past maxtime = $(maxtime)"))
			end

		#	Map each day to the latest log row at or before it
			row_for_day = Vector{Int}(undef, maxtime + 1)
			k = 1
			for day in 0:maxtime
				while k < nrow(ilog) && times[k + 1] <= day
					k += 1
				end
				row_for_day[day + 1] = k
			end

		#	Build the daily table and set the day column
			daily = ilog[row_for_day, :]
			daily.time = Float64.(0:maxtime)

		#	Assembling Result
			return daily
	end
	@doc raw"""
	**Description**
	Expand a `sirdif` infection log to exactly one row per day from 0 to `maxtime`, filling the days after extinction.

	**Usage**
	`expand_daily_log(ilog, maxtime)`

	**Arguments**
	- `ilog::DataFrame`: The `infection_log` returned by `sirdif`.
	- `maxtime::Int`: The simulation horizon in days.

	**Details**
	The simulator writes one row per day until the last infection clears, then a single row at `maxtime`. Each missing day copies the latest row at or before it. After extinction no further infection is possible, so the infected count and cumulative incidence are constant, and carrying them forward is exact. The simulator log itself is left unchanged, which keeps the SAS comparison intact.

	**Value**
	A `DataFrame` with `maxtime + 1` rows, days `0:maxtime`, and the same columns as `ilog`.

	**See Also**
	`trajectory_summaries`, `summarize_run`
	""" expand_daily_log

# ============================================================
#	TRAJECTORY SUMMARIES
# ============================================================

#	Helper Function for trajectory_summaries: first day a share of growth is reached
	function growth_crossing_day(cum_t::Vector{Float64}, share::Float64)::Float64
		"""
		Args:
			cum_t::Vector{Float64}: proportion ever infected on days 0, 1, ..., T
			share::Float64: share of the growth beyond day 0, in (0, 1]
		Returns:
			Float64: first day d with cum_t[d] - cum_t[0] ≥ share × (cum_t[T] - cum_t[0]),
					 or NaN when there is no growth (a non-starter)
		Notes:
			Growth is measured beyond day 0, so the seeds do not count,
			matching Option A in the Richards fit. Because cum_t is
			non-decreasing and the final value meets every share, a day is
			always found when there is growth.
		"""

		#	No growth beyond the seeds
			growth = cum_t[end] - cum_t[1]
			if growth <= 0.0
				return NaN
			end

		#	First day at or above the share of growth
			idx = findfirst(c -> c - cum_t[1] >= share * growth, cum_t)
			return Float64(idx - 1)
	end

#	Trajectory Summaries from a Daily Series
	function trajectory_summaries(infected_t::Vector{Float64}, cum_t::Vector{Float64};
								  surge_time::Int=21)
		"""
		Args:
			infected_t::Vector{Float64}: proportion currently infected on days 0, 1, ..., T
			cum_t::Vector{Float64}: proportion ever infected on days 0, 1, ..., T
			surge_time::Int: day of the early-infection measurement (default = 21)
		Returns:
			NamedTuple: peak_infected, peak_day, infected_at_surge, max_slope,
						final_size, duration, extinguished, t_10, t_half, t_90
		Notes:
			peak_day is the first day of the peak. infected_at_surge is NaN if
			surge_time lies beyond the series. max_slope is the largest
			one-day increase in the proportion infected. duration is the first
			day after day 0 on which no one is infected; if infections persist
			through the horizon, duration is the last day and extinguished is
			false (the duration is censored).
			t_10, t_half, and t_90 are the first days on which the run reaches
			10%, 50%, and 90% of its growth beyond the seeds (26 October 2026).
			They come from the data, not the Richards fit. For runs whose fit
			reaches a = -1, the fitted t_i marks the takeoff rather than the
			peak, and these days keep a consistent meaning in every run.
			t_90 - t_half measures the length of the tail. All three are NaN
			for a non-starter.
		"""

		#	Validation
			if isempty(infected_t) || length(infected_t) != length(cum_t)
				throw(ArgumentError("series must be non-empty and of equal length"))
			end

		#	Peak, surge, and slope
			peak_idx = argmax(infected_t)
			peak_infected = infected_t[peak_idx]
			infected_at_surge = surge_time + 1 <= length(infected_t) ? infected_t[surge_time + 1] : NaN
			max_slope = length(infected_t) > 1 ? maximum(diff(infected_t)) : 0.0

		#	Duration: first day after day 0 with no one infected
			end_idx = findfirst(==(0.0), infected_t[2:end])
			extinguished = end_idx !== nothing
			duration = extinguished ? end_idx : length(infected_t) - 1

		#	Timing of growth: 10%, 50%, and 90% of the growth beyond the seeds
			t_10 = growth_crossing_day(cum_t, 0.1)
			t_half = growth_crossing_day(cum_t, 0.5)
			t_90 = growth_crossing_day(cum_t, 0.9)

		#	Assembling Result
			return (peak_infected=peak_infected, peak_day=peak_idx - 1, infected_at_surge=infected_at_surge,
					max_slope=max_slope, final_size=cum_t[end], duration=duration, extinguished=extinguished,
					t_10=t_10, t_half=t_half, t_90=t_90)
	end
	@doc raw"""
	**Description**
	Compute the trajectory-level planning summaries for one run from its daily series.

	**Usage**
	`trajectory_summaries(infected_t, cum_t; surge_time=21)`

	**Arguments**
	- `infected_t::Vector{Float64}`: Proportion currently infected on days 0 to T.
	- `cum_t::Vector{Float64}`: Proportion ever infected on days 0 to T.
	- `surge_time::Int`: Day of the early-infection measurement (default 21).

	**Details**
	The three timing days are taken from the data rather than the Richards fit. When the fit reaches its lower shape limit, ``a = -1``, the fitted ``t_i`` falls on the day of takeoff rather than the peak of incidence, because the monomolecular curve jumps to its largest incidence at its start. In simulated village epidemics this happens for runs with an abrupt takeoff, sometimes after a delay, followed by a long tail with secondary waves. ``t_{\text{half}}`` keeps the same meaning in every run, and ``t_{90} - t_{\text{half}}`` measures the length of the tail.

	**Value**
	A `NamedTuple`:
	- `peak_infected`: Largest proportion simultaneously infected.
	- `peak_day`: First day on which that peak occurs.
	- `infected_at_surge`: Proportion infected on day `surge_time` (`NaN` beyond the series).
	- `max_slope`: Largest one-day increase in the proportion infected.
	- `final_size`: Proportion ever infected on the last day, C(T).
	- `duration`: First day after day 0 with no one infected, or T if infection persists.
	- `extinguished`: `false` when infection persists through day T, so `duration` is censored.
	- `t_10`, `t_half`, `t_90`: First days on which the run reaches 10%, 50%, and 90% of its growth beyond the seeds (`NaN` for a non-starter).

	**See Also**
	`expand_daily_log`, `summarize_run`
	""" trajectory_summaries

# ============================================================
#	RUN SUMMARY
# ============================================================

#	Summarize One Simulated Run
	function summarize_run(ilog::DataFrame, n_nodes::Int, maxtime::Int, region::TargetRegion;
						   rec_scale::Union{Float64, Nothing}=nothing,
						   rec_shape::Union{Float64, Nothing}=nothing,
						   extinct_policy::Symbol=:exclude,
						   extinction_threshold::Float64=0.02)
		"""
		Args:
			ilog::DataFrame: infection_log returned by sirdif
			n_nodes::Int: number of nodes in the network
			maxtime::Int: simulation horizon in days
			region::TargetRegion: acceptable outcome region
			rec_scale, rec_shape: Weibull recovery parameters of the run
								  (needed for peak_prevalence and P_max_thresh)
			extinct_policy::Symbol: :exclude, :success, or :failure (default = :exclude)
			extinction_threshold::Float64: passed to fit_richards (default = 0.02)
		Returns:
			NamedTuple: Richards fit and its residual error, peak quantities,
						trajectory summaries, run classification, and target evaluation
		Notes:
			Pipeline: expand to one row per day, read C0 from day 0, fit the
			Richards curve to cumulative incidence (Option A), compute the
			trajectory summaries, and evaluate the target region.
			C0 is the day-0 cumulative incidence, which sirdif records as
			(number of distinct seeds) / N, so it is exact.
			The infected proportion is n_infected / n_nodes. The log's
			prop_currently_infected is 0.0 on day 0 by SAS convention even
			though the seeds are infected, so it is not used.
			t_i is the fitted inflection time, the day of peak incidence. With
			K, r_growth, and a it reproduces the fitted curve, and fit_rmse
			reports how closely that curve matches the run.
			a lies in [-1, 10]: 1 is logistic, 0 Gompertz, and -1
			monomolecular. frac_at_inflection, (1 + a)^(-1/a), is the share of
			K reached at t_i; a value below 1/e marks a front-loaded epidemic
			with an early peak and a long tail.
		"""

		#	Validation
			if n_nodes < 1
				throw(ArgumentError("n_nodes must be ≥ 1"))
			end
			if !("prop_cum_infected" in names(ilog))
				throw(ArgumentError("infection log lacks prop_cum_infected; use the revised sirdif"))
			end

		#	Expand to one row per day
			daily = expand_daily_log(ilog, maxtime)
			t = Float64.(daily.time)
			cum_t = Float64.(daily.prop_cum_infected)
			infected_t = Float64.(daily.n_infected) ./ n_nodes
			C0 = cum_t[1]

		#	Fit the Richards curve beyond the seeds
			fit = fit_richards(t, cum_t; C0=C0, extinction_threshold=extinction_threshold)

		#	Trajectory summaries
			traj = trajectory_summaries(infected_t, cum_t; surge_time=region.surge_time)

		#	Analytic peak prevalence when the fit and recovery parameters allow it
			if fit.fit_ok && rec_scale !== nothing && rec_shape !== nothing
				p_max = peak_prevalence(fit.K, fit.r_growth, fit.a, rec_scale, rec_shape)
			else
				p_max = NaN
			end

		#	Classify the run and evaluate the target
			target = evaluate_target(region, infected_t, fit.K, fit.r_growth, fit.t_i, fit.a;
									 C0=C0, rec_scale=rec_scale, rec_shape=rec_shape,
									 status=fit.status, extinct_policy=extinct_policy,
									 final_cum=cum_t[end])

		#	Assembling Result
			return (K=fit.K, r_growth=fit.r_growth, t_i=fit.t_i, a=fit.a,
					frac_at_inflection=fit.frac_at_inflection, fit_rmse=fit.fit_rmse,
					peak_incidence=fit.peak_incidence, peak_prevalence=p_max,
					peak_infected=traj.peak_infected, peak_day=traj.peak_day,
					infected_at_surge=traj.infected_at_surge, max_slope=traj.max_slope,
					final_size=traj.final_size, duration=traj.duration, extinguished=traj.extinguished,
					t_10=traj.t_10, t_half=traj.t_half, t_90=traj.t_90, C0=C0, fit_ok=fit.fit_ok, extinct=fit.extinct, viable=fit.viable,
					status=fit.status, class=target.class, in_target=target.in_target,
					peak_ok=target.peak, surge_ok=target.surge, slope_ok=target.slope,
					prevalence_ok=target.prevalence, final_size_ok=target.final_size,
					n_days=nrow(daily))
	end
	@doc raw"""
	**Description**
	Turn one `sirdif` infection log into the summary record used for inference: the Richards fit and its residual error, peak quantities, trajectory summaries, run classification, and target evaluation.

	**Usage**
	```julia
	summarize_run(ilog, n_nodes, maxtime, region; rec_scale=nothing, rec_shape=nothing,
				  extinct_policy=:exclude, extinction_threshold=0.02)
	```

	**Arguments**
	- `ilog::DataFrame`: The `infection_log` returned by the revised `sirdif` (with `prop_cum_infected`).
	- `n_nodes::Int`: Number of nodes in the network.
	- `maxtime::Int`: Simulation horizon in days.
	- `region::TargetRegion`: The acceptable outcome region.
	- `rec_scale, rec_shape`: The run's Weibull recovery parameters; required when `region.P_max_thresh` is set, and used for `peak_prevalence`.
	- `extinct_policy::Symbol`: Treatment of extinct runs (default `:exclude`).
	- `extinction_threshold::Float64`: Growth beyond the seeds below which a run is extinct (default `0.02`).

	**Details**
	1. `expand_daily_log` fills the days after extinction.
	2. ``C_0`` is read from day 0, which `sirdif` records as the number of distinct seeds divided by ``N``.
	3. `fit_richards` fits ``C(t) - C_0`` and classifies the run (`:ok`, `:nonstarter`, `:extinct`, or a failure status). The fitted `t_i` is the inflection time, the day of peak incidence; with `K`, `r_growth`, and `a` it reproduces the fitted curve. The shape `a` lies in ``[-1, 10]``, and `frac_at_inflection` gives the share of `K` reached at `t_i`.
	4. `trajectory_summaries` computes peak, early-infection size, slope, final size, duration, and the days on which 10%, 50%, and 90% of the growth beyond the seeds is reached, from ``I(t) = n_{\text{infected}}(t)/N``.
	5. `evaluate_target` assigns the class and `in_target` under `extinct_policy`.

	Every analysis is conditioned on at least one transmission, so a non-starter has `viable = false` and `in_target = nothing`.

	**Value**
	A `NamedTuple` with the fitted parameters (`K`, `r_growth`, `t_i`, `a`), the shape summary (`frac_at_inflection`), the fit's root-mean-square residual (`fit_rmse`), peak quantities (`peak_incidence`, `peak_prevalence`), trajectory summaries (`peak_infected`, `peak_day`, `infected_at_surge`, `max_slope`, `final_size`, `duration`, `extinguished`), the timing of growth (`t_10`, `t_half`, `t_90`), `C0`, the flags (`fit_ok`, `extinct`, `viable`, `status`, `class`), the target result (`in_target` and each constraint), and `n_days`.

	**See Also**
	`expand_daily_log`, `trajectory_summaries`, `fit_richards`, `evaluate_target`
	""" summarize_run

#   Exporting Objects
	export expand_daily_log,
		   trajectory_summaries,
		   summarize_run

end # module run_summary

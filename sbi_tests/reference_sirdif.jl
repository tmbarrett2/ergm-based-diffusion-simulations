#Reference Copy of the Validated sirdif
#Frozen verbatim from diffusion_sim.jl as of 29 September 2026

"""
Frozen copy of the validated sirdif() used as the regression baseline
for UT-1.7. Do not edit: any change breaks the meaning of the test.
The nested remove_id! definition is kept verbatim for that reason.
"""
module ReferenceSirdif

using DataFrames
using Random

#	SIR Diffusion Simulation with Social Feedback
	function sirdif(
		alst::Union{Matrix{Int}, Vector{Vector{Int}}},
		vlst::Union{Matrix{Float64}, Vector{Vector{Float64}}},
		infectedp::Vector{Int}, inf_r::Float64, rec_t::Int,
		maxtime::Int, p_symp::Float64, b_int::Float64,
		b_close::Float64, b_cxn_peers::Float64, b_cxn_total::Float64,
		b_cxn_symp::Float64, b_cls_x_smp::Float64;
		transmission_method::Symbol = :weighted,
		immunity_duration::Union{Int,Nothing} = nothing
	)
		"""
		Args:
			alst: Adjacency list (Matrix or Vector{Vector})
			vlst: Edge weights aligned to alst
			infectedp: Initial infected node IDs
			inf_r: Base transmission probability
			rec_t: Recovery time in days
			maxtime: Maximum simulation days
			p_symp: Probability infected is symptomatic
			b_int: Baseline interaction coefficient
			b_close: Edge weight coefficient
			b_cxn_peers: Peer infection coefficient
			b_cxn_total: Global infection coefficient
			b_cxn_symp: Symptomatic coefficient
			b_cls_x_smp: Peers × symptomatic interaction
			transmission_method: :weighted or :simple (default :weighted)
			immunity_duration: Days of immunity after recovery (default nothing)
		Returns:
			Dict with infection_log, total_time, final_state
		Notes:
			Exact SAS behavior: FIFO ordering, cumulative peer counts, and single-row padding at extinction.
		"""

		#	Cache network size
			if alst isa Matrix{Int}
				n_nodes = size(alst, 1)
			else
				n_nodes = length(alst)
			end

		#	Normalize inputs to vector-of-vectors; collect node ids
			if alst isa Matrix{Int}
				alst_vec   = Vector{Vector{Int}}(undef, n_nodes)
				vlst_vec   = Vector{Vector{Float64}}(undef, n_nodes)
				unique_ids = Vector{Int}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					ego_id        = alst[i, 1]
					unique_ids[i] = ego_id
					row_ids       = alst[i, 2:end]
					row_wgts      = vlst[i, 2:end]
					nz_mask       = row_ids .> 0
					neighbors     = row_ids[nz_mask]
					weights       = row_wgts[nz_mask]
					alst_vec[i]   = [ego_id; neighbors]
					vlst_vec[i]   = [Float64(ego_id); weights]
				end
			else
				alst_vec   = Vector{Vector{Int}}(undef, n_nodes)
				vlst_vec   = Vector{Vector{Float64}}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					alst_vec[i] = alst[i]
					vlst_vec[i] = vlst[i]
				end
				unique_ids = Vector{Int}(undef, n_nodes)
				@inbounds for i in 1:n_nodes
					unique_ids[i] = alst_vec[i][1]
				end
			end

		#	Map node id → row index
			id_to_idx = Dict{Int,Int}(unique_ids[i] => i for i in 1:n_nodes)

		#	Start wall clock
			start_time = time()

		#	State matrix: [id, I, S, R, t_rec, nbrsinf, infection_order]
			state = zeros(Int, n_nodes, 7)
			@inbounds begin
				state[:, 1] .= unique_ids
				state[:, 3] .= 1
			end

		#	Column indices
			ID_COL               = 1
			INFECTED_COL         = 2
			SUSCEPTIBLE_COL      = 3
			RECOVERED_COL        = 4
			TIME_TO_RECOVERY_COL = 5
			NBRSINF_COL          = 6
			INFECTION_ORDER_COL  = 7

		#	Initialize infection order counter
			infection_counter = 0

		#	Seed infections
			@inbounds for inf_id in infectedp
				idx = id_to_idx[inf_id]
				state[idx, INFECTED_COL]         = 1
				state[idx, SUSCEPTIBLE_COL]      = 0
				state[idx, TIME_TO_RECOVERY_COL] = rec_t
				infection_counter += 1
				state[idx, INFECTION_ORDER_COL]  = infection_counter
			end

		#	Susceptible adjacency (preserve original column order)
			s_alst_vec = Vector{Vector{Int}}(undef, n_nodes)
			s_vlst_map = Vector{Dict{Int,Float64}}(undef, n_nodes)
			@inbounds for i in 1:n_nodes
				neighbors   = alst_vec[i][2:end]
				weights     = vlst_vec[i][2:end]
				s_alst_vec[i] = copy(neighbors)
				s_vlst_map[i] = Dict(zip(neighbors, weights))
			end

		#	In-place vector filter helper
			@inline function remove_id!(v::Vector{Int}, id::Int)
				w = 1
				@inbounds for i in 1:length(v)
					x = v[i]
					if x != id
						v[w] = x
						w += 1
					end
				end
				if w <= length(v)
					resize!(v, w - 1)
				end
				return nothing
			end

		#	Remove initially infected from all susceptible lists
			@inbounds for inf_id in infectedp
				for i in 1:n_nodes
					remove_id!(s_alst_vec[i], inf_id)
				end
			end

		#	Persistent peer-infected counts (cumulative)
			nbrsinf = zeros(Int, n_nodes)
			@inbounds for i in 1:n_nodes
				if state[i, INFECTED_COL] == 1
					for alter_id in s_alst_vec[i]
						alter_idx = id_to_idx[alter_id]
						nbrsinf[alter_idx] += 1
					end
				end
			end

		#	Time series buffer (preallocated; may finish early)
			timesum = Matrix{Float64}(undef, maxtime + 1, 7)
			wrow    = 1
			n_initial = length(infectedp)
			pinf      = n_initial / n_nodes
			timesum[wrow, :] = Float64[0.0, n_initial, pinf, 0.0, 0.0, 0.0, 0.0]

		#	Per-day buffers
			infected_idx_buf   = Vector{Int}(undef, n_nodes)
			infected_order_buf = Vector{Int}(undef, n_nodes)
			infected_perm_buf  = Vector{Int}(undef, 0)
			ordered_idx_buf    = Vector{Int}(undef, n_nodes)

		#	Main loop over days
			@inbounds for t in 1:maxtime

				#	Collect infected at start-of-day
					ninf = 0
					for i in 1:n_nodes
						if state[i, INFECTED_COL] == 1
							ninf += 1
							infected_idx_buf[ninf]   = i
							infected_order_buf[ninf] = state[i, INFECTION_ORDER_COL]
						end
					end

				#	Extinction padding (SAS: single row at maxtime)
					if ninf == 0
						NR_now = 0
						for i in 1:n_nodes
							NR_now += state[i, RECOVERED_COL]
						end
						NI_now = 0
						prop_ever_now = NR_now == 0 ? 0.0 : NR_now / n_nodes
						prop_cur_now  = 0.0
						prop_rec_now  = NR_now == 0 ? 0.0 : 1.0
						wrow += 1
						timesum[wrow, :] = Float64[maxtime, NI_now, prop_ever_now, prop_cur_now, NR_now, prop_rec_now, NR_now]
						break
					end

				#	Order infected FIFO by infection_order
					resize!(infected_perm_buf, ninf)
					for i in 1:ninf
						infected_perm_buf[i] = i
					end
					sort!(infected_perm_buf; by = i -> infected_order_buf[i])
					for ii in 1:ninf
						ordered_idx_buf[ii] = infected_idx_buf[infected_perm_buf[ii]]
					end

				#	Process each infected ego
					for p in 1:ninf
						ego_idx = ordered_idx_buf[p]

						#	Neighbors in original column order
							susceptible_neighbor_ids = s_alst_vec[ego_idx]

						#	Ego symptomatic draw (per-day)
							issympt = (rand() < p_symp) ? 1 : 0

						#	Scan neighbors
							for alter_id in susceptible_neighbor_ids
								alter_idx = id_to_idx[alter_id]
								edgwgt    = s_vlst_map[ego_idx][alter_id]

								#	Current cumulative peers-infected
									peersinf = nbrsinf[alter_idx]

								#	Activation logit
									lwact = b_int +
									        (issympt * b_cxn_symp) +
									        (b_close * edgwgt) +
									        (b_cxn_peers * peersinf) +
									        (b_cxn_total * (ninf / (n_nodes / 3))) +
									        (b_cls_x_smp * peersinf * issympt)
									prob_act = exp(lwact) / (1 + exp(lwact))

								#	Activation and transmission
									if rand() < prob_act
										transprob = transmission_method === :weighted ? (1 - (1 - inf_r)^edgwgt) : inf_r
										if rand() < transprob
											#	Update state
												state[alter_idx, INFECTED_COL]         = 1
												state[alter_idx, SUSCEPTIBLE_COL]      = 0
												state[alter_idx, TIME_TO_RECOVERY_COL] = rec_t
												infection_counter += 1
												state[alter_idx, INFECTION_ORDER_COL]  = infection_counter

											#	Remove alter from all susceptible lists
												for k in 1:n_nodes
													remove_id!(s_alst_vec[k], alter_id)
												end

											#	Increment peers-infected for ego’s original neighbors
												for nbr_id in @view alst_vec[ego_idx][2:end]
													nbr_idx = id_to_idx[nbr_id]
													nbrsinf[nbr_idx] += 1
												end
										end
									end
							end

						#	Recovery countdown
							state[ego_idx, TIME_TO_RECOVERY_COL] -= 1
					end

				#	Move recovered
					for i in 1:n_nodes
						if state[i, TIME_TO_RECOVERY_COL] == 0 && state[i, INFECTED_COL] == 1
							state[i, INFECTED_COL]         = 0
							state[i, RECOVERED_COL]        = 1
							state[i, TIME_TO_RECOVERY_COL] = 0
							state[i, INFECTION_ORDER_COL]  = 0
							if immunity_duration !== nothing
								state[i, TIME_TO_RECOVERY_COL] = -immunity_duration
							end
						end
					end

				#	SIRS waning immunity (optional)
					if immunity_duration !== nothing
						for i in 1:n_nodes
							if state[i, TIME_TO_RECOVERY_COL] < 0 && state[i, RECOVERED_COL] == 1
								state[i, TIME_TO_RECOVERY_COL] += 1
							end
						end
						for i in 1:n_nodes
							if state[i, TIME_TO_RECOVERY_COL] == 0 && state[i, RECOVERED_COL] == 1
								state[i, RECOVERED_COL]   = 0
								state[i, SUSCEPTIBLE_COL] = 1
								for j in 1:n_nodes
									if i != j && (unique_ids[i] in alst_vec[j][2:end])
										push!(s_alst_vec[j], unique_ids[i])	# Only in SIRS mode
									end
								end
							end
						end
					end

				#	State consistency check
					for i in 1:n_nodes
						sum_row = state[i, INFECTED_COL] + state[i, SUSCEPTIBLE_COL] + state[i, RECOVERED_COL]
						@assert sum_row == 1 "Invalid state for node $(state[i, ID_COL]): I=$(state[i,2]) S=$(state[i,3]) R=$(state[i,4])"
					end

				#	Record daily metrics
					NI = 0
					NR = 0
					for i in 1:n_nodes
						NI += state[i, INFECTED_COL]
						NR += state[i, RECOVERED_COL]
					end
					prop_ever = (NI + NR) / n_nodes
					prop_cur  = NI / n_nodes
					prop_rec  = (NI + NR) > 0 ? NR / (NI + NR) : 0.0

					wrow += 1
					timesum[wrow, :] = Float64[t, NI, prop_ever, prop_cur, NR, prop_rec, NR]
			end

		#	Stop clock
			total_time = time() - start_time

		#	Trim time series
			inflog = timesum[1:wrow, 1:6]
			inflog_df = DataFrame(time = inflog[:, 1], n_infected = inflog[:, 2], prop_ever_infected = inflog[:, 3], 
								  prop_currently_infected = inflog[:, 4], n_recovered = inflog[:, 5], prop_recovered = inflog[:, 6])

		#	Return results
			return Dict{String,Any}(
				"infection_log" => inflog_df,
				"total_time"    => total_time,
				"final_state"   => state,
			)
	end

end # module ReferenceSirdif

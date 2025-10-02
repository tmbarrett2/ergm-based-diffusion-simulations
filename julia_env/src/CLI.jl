module CLI

using ArgParse
using CSV
using JSON3
using Random
using DataFrames
using ProgressMeter
using ..diffusion_sim: sim_prep, sirdif

#	Helper Function for CLI: write results (CSV + JSON metadata)
	function _write_results(payload, outpath::AbstractString = "-")
		"""
		Args:
			payload: Dict with infection_log, total_time, final_state/final_states
			outpath: Base path for output files (- for stdout)
		Returns:
			Nothing
		Notes:
			Handles multi-network metadata if network_id column exists
		"""
		
		#	Extract components
			infection_log = payload["infection_log"]
			
		#	Check if multi-network (has network_id column)
			is_multi = "network_id" in names(infection_log)
			
		#	Build metadata
			if is_multi
				#	Group by network and get final state for each
					networks = unique(infection_log.network_id)
					network_summaries = Dict[]
					
					for net_id in networks
						net_data = filter(row -> row.network_id == net_id, infection_log)
						final_row = net_data[end, :]
						
						push!(network_summaries, Dict(
							"network_id" => net_id,
							"final_n_infected" => final_row.n_infected,
							"final_n_recovered" => final_row.n_recovered,
							"final_prop_ever_infected" => final_row.prop_ever_infected,
							"max_time_reached" => maximum(net_data.time)
						))
					end
					
					metadata = Dict(
						"total_networks" => length(networks),
						"network_summaries" => network_summaries,
						"total_time" => payload["total_time"]
					)
			else
				#	Single network metadata
					metadata = Dict(
						"total_time" => payload["total_time"],
						"final_state_summary" => Dict(
							"n_susceptible" => sum(payload["final_state"][:, 3]),
							"n_infected"    => sum(payload["final_state"][:, 2]),
							"n_recovered"   => sum(payload["final_state"][:, 4])
						)
					)
			end

		#	Write based on output path
			if outpath == "-"
				CSV.write(stdout, infection_log)
			else
				csv_path  = replace(outpath, r"\.json$" => "") * "_infection_log.csv"
				CSV.write(csv_path, infection_log)

				json_path = replace(outpath, r"\.json$" => "") * "_metadata.json"
				open(json_path, "w") do io
					write(io, JSON3.write(metadata))
				end

				println("Results written to:")
				println("  Infection log: $csv_path")
				println("  Metadata: $json_path")
			end

		#	Return
			return nothing
	end

#	Single Wrapper: sim_prep → sirdif pipeline
	function run_simulation(
		graphml::AbstractString;
		sas_transformation::Bool = false,
		infectedp::Vector{Int},
		inf_r::Float64, rec_t::Int, maxtime::Int, p_symp::Float64,
		b_int::Float64, b_close::Float64, b_cxn_peers::Float64, b_cxn_total::Float64,
		b_cxn_symp::Float64, b_cls_x_smp::Float64,
		transmission_method::Symbol = :weighted,
		immunity_duration::Union{Int,Nothing} = nothing,
		seed::Union{Int,Nothing} = nothing)
		"""
		Args:
			graphml: Path to GraphML file
			sas_transformation: Apply SAS preprocessing
			infectedp: Initial infected node IDs
			inf_r: Base transmission probability
			rec_t: Recovery time
			maxtime: Maximum days
			p_symp: Symptomatic probability
			b_*: Logit coefficients
			transmission_method: :weighted or :simple
			immunity_duration: SIRS immunity days
			seed: RNG seed
		Returns:
			Dict with infection_log, total_time, final_state
		"""

		#	Validate inputs
			if !isfile(graphml)
				throw(ArgumentError("GraphML path not found: $graphml"))
			end
			if !(transmission_method === :weighted || transmission_method === :simple)
				throw(ArgumentError("transmission_method must be :weighted or :simple"))
			end

		#	Seed (optional)
			if seed !== nothing
				Random.seed!(seed)
			end

		#	Prepare network
			data = sim_prep(graphml; sas_transformation = sas_transformation)
			alst = data.alst
			vlst = data.vlst

		#	Run simulation
			res = sirdif(
				alst, vlst, infectedp, inf_r, rec_t, maxtime, p_symp,
				b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp;
				transmission_method = transmission_method,
				immunity_duration = immunity_duration
			)

		#	Return Dict
			return res
	end

#	Sweep Wrapper: parameter grid sweep
	function run_sweep(
		graphml::AbstractString;
		sas_transformation::Bool = false,
		infectedp::Vector{Int},
		inf_r::Float64, rec_t::Int, maxtime::Int, p_symp::Float64,
		b_int::Float64, b_close::Float64,
		b_cxn_peers::Vector{Float64}, b_cxn_total::Vector{Float64}, b_cxn_symp::Vector{Float64},
		b_cls_x_smp::Float64,
		num_iter::Int = 1,
		transmission_method::Symbol = :weighted,
		immunity_duration::Union{Int,Nothing} = nothing,
		seed::Union{Int,Nothing} = nothing)
		"""
		Args:
			graphml: Path to GraphML file
			b_cxn_peers/total/symp: Vectors for parameter sweep
			num_iter: Iterations per parameter combo
			[other args same as run_simulation]
		Returns:
			Dict with stacked infection_log DataFrame
		"""

		#	Optional global seed
			if seed !== nothing
				Random.seed!(seed)
			end

		#	Load network once
			data = sim_prep(graphml; sas_transformation = sas_transformation)
			alst = data.alst
			vlst = data.vlst

		#	Prepare output + progress
			results     = DataFrame()
			total_runs  = length(b_cxn_peers) * length(b_cxn_total) * length(b_cxn_symp) * num_iter
			progressbar = Progress(total_runs; desc = "Parameter sweep", dt = 0.2)

		#	Triple loop over coefficients + iterations
			local last_res = nothing
			for bp in b_cxn_peers
				for bg in b_cxn_total
					for bs in b_cxn_symp
						for it in 1:num_iter
							last_res = sirdif(
								alst, vlst, infectedp, inf_r, rec_t, maxtime, p_symp,
								b_int, b_close, bp, bg, bs, b_cls_x_smp;
								transmission_method = transmission_method,
								immunity_duration = immunity_duration
							)

							ilog = last_res["infection_log"]
							n = nrow(ilog)
							ilog[!, :b_cxn_peers] = fill(bp, n)
							ilog[!, :b_cxn_total] = fill(bg, n)
							ilog[!, :b_cxn_symp]  = fill(bs, n)
							ilog[!, :iter]        = fill(it, n)

							append!(results, ilog)
							next!(progressbar)
						end
					end
				end
			end

		#	Return stacked results
			return Dict(
				"infection_log" => results,
				"total_time"    => maxtime,
				"final_state"   => last_res === nothing ? zeros(Int, 0, 0) : last_res["final_state"]
			)
	end

#	Multi-Network Wrapper: process multiple GraphML files
	function run_multi(
		graphml_files::Vector{String};
		sas_transformation::Bool = false,
		infectedp::Vector{Int},
		inf_r::Float64, rec_t::Int, maxtime::Int, p_symp::Float64,
		b_int::Float64, b_close::Float64, b_cxn_peers::Float64, b_cxn_total::Float64,
		b_cxn_symp::Float64, b_cls_x_smp::Float64,
		transmission_method::Symbol = :weighted,
		immunity_duration::Union{Int,Nothing} = nothing,
		seed::Union{Int,Nothing} = nothing)
		"""
		Args:
			graphml_files: Vector of GraphML file paths
			[other args same as run_simulation]
		Returns:
			Dict with combined infection_log DataFrame and final_states
		Notes:
			Adds network_id column extracted from filename
			Stores final_state for each network
		"""
		
		#	Optional global seed
			if seed !== nothing
				Random.seed!(seed)
			end

		#	Prepare output
			results = DataFrame()
			final_states = Dict{String, Matrix{Int}}()
			progressbar = Progress(length(graphml_files); desc = "Processing networks", dt = 0.2)

		#	Process each network
			for graphml in graphml_files
				#	Extract network ID from path
					network_id = replace(basename(graphml), r"\.graphml$" => "")
				
				#	Load and run
					data = sim_prep(graphml; sas_transformation = sas_transformation)
					res = sirdif(
						data.alst, data.vlst, infectedp, inf_r, rec_t, maxtime, p_symp,
						b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp;
						transmission_method = transmission_method,
						immunity_duration = immunity_duration
					)

				#	Store final state for this network
					final_states[network_id] = res["final_state"]

				#	Add network identifier
					ilog = res["infection_log"]
					insertcols!(ilog, 1, :network_id => network_id)

				#	Append
					append!(results, ilog)
					next!(progressbar)
			end

		#	Return combined results with all final states
			return Dict(
				"infection_log" => results,
				"total_time"    => maxtime,
				"final_states"  => final_states,
				"final_state"   => zeros(Int, 0, 0)  # Empty for compatibility
			)
	end

#	CLI Main Function (modes: single | sweep | multi)
	function cli_main(args::Vector{String} = ARGS)
		"""
		Args:
			args: Command line arguments
		Returns:
			Nothing
		Notes:
			Three modes: single run, parameter sweep, multi-network
		"""
		
		#	No args: print usage
			if isempty(args)
				println("Usage:")
				println("  single: diffusion_sim.CLI single --graphml <file> ... --out <path>")
				println("  sweep:  diffusion_sim.CLI sweep  --graphml <file> ... --out <path>")
				println("  multi:  diffusion_sim.CLI multi  --graphml <file1,file2> ... --out <path>")
				return nothing
			end

		#	Mode select
			mode = args[1]

		#	Schema: SINGLE
			if mode == "single"
				s = ArgParseSettings()
				@add_arg_table! s begin
					"--graphml"
						arg_type = String
						required = true
					"--sas-transformation"
						action = :store_true
					"--infected"
						arg_type = Int
						nargs = '+'
						required = true
					"--inf-r"
						arg_type = Float64
						required = true
					"--rec-t"
						arg_type = Int
						required = true
					"--max-time"
						arg_type = Int
						required = true
					"--p-symp"
						arg_type = Float64
						required = true
					"--b-int"
						arg_type = Float64
						required = true
					"--b-close"
						arg_type = Float64
						required = true
					"--b-cxn-peers"
						arg_type = Float64
						required = true
					"--b-cxn-total"
						arg_type = Float64
						required = true
					"--b-cxn-symp"
						arg_type = Float64
						required = true
					"--b-cls-x-smp"
						arg_type = Float64
						required = true
					"--transmission"
						arg_type = String
						default = "weighted"
					"--immunity-duration"
						arg_type = Int
					"--seed"
						arg_type = Int
					"--out"
						arg_type = String
						default = "-"
				end

				p = parse_args(args[2:end], s)
				if p === nothing; return nothing; end

				res = run_simulation(
					p["graphml"];
					sas_transformation = get(p, "sas-transformation", false),
					infectedp = Vector{Int}(p["infected"]),
					inf_r = p["inf-r"],
					rec_t = p["rec-t"],
					maxtime = p["max-time"],
					p_symp = p["p-symp"],
					b_int = p["b-int"],
					b_close = p["b-close"],
					b_cxn_peers = p["b-cxn-peers"],
					b_cxn_total = p["b-cxn-total"],
					b_cxn_symp = p["b-cxn-symp"],
					b_cls_x_smp = p["b-cls-x-smp"],
					transmission_method = p["transmission"] == "simple" ? :simple : :weighted,
					immunity_duration = get(p, "immunity-duration", nothing),
					seed = get(p, "seed", nothing)
				)
				_write_results(res, p["out"])
				return nothing
			end

		#	Schema: SWEEP
			if mode == "sweep"
				s = ArgParseSettings()
				@add_arg_table! s begin
					"--graphml"
						arg_type = String
						required = true
					"--sas-transformation"
						action = :store_true
					"--infected"
						arg_type = Int
						nargs = '+'
						required = true
					"--inf-r"
						arg_type = Float64
						required = true
					"--rec-t"
						arg_type = Int
						required = true
					"--max-time"
						arg_type = Int
						required = true
					"--p-symp"
						arg_type = Float64
						required = true
					"--b-int"
						arg_type = Float64
						required = true
					"--b-close"
						arg_type = Float64
						required = true
					"--b-cxn-peers"
						arg_type = String
						required = true
					"--b-cxn-total"
						arg_type = String
						required = true
					"--b-cxn-symp"
						arg_type = String
						required = true
					"--b-cls-x-smp"
						arg_type = Float64
						required = true
					"--num-iter"
						arg_type = Int
						default = 1
					"--transmission"
						arg_type = String
						default = "weighted"
					"--immunity-duration"
						arg_type = Int
					"--seed"
						arg_type = Int
					"--out"
						arg_type = String
						default = "-"
				end

				p = parse_args(args[2:end], s)
				if p === nothing; return nothing; end

				bp = parse.(Float64, strip.(split(p["b-cxn-peers"], ',')))
				bg = parse.(Float64, strip.(split(p["b-cxn-total"], ',')))
				bs = parse.(Float64, strip.(split(p["b-cxn-symp"], ',')))

				res = run_sweep(
					p["graphml"];
					sas_transformation = get(p, "sas-transformation", false),
					infectedp = Vector{Int}(p["infected"]),
					inf_r = p["inf-r"],
					rec_t = p["rec-t"],
					maxtime = p["max-time"],
					p_symp = p["p-symp"],
					b_int = p["b-int"],
					b_close = p["b-close"],
					b_cxn_peers = bp,
					b_cxn_total = bg,
					b_cxn_symp = bs,
					b_cls_x_smp = p["b-cls-x-smp"],
					num_iter = p["num-iter"],
					transmission_method = p["transmission"] == "simple" ? :simple : :weighted,
					immunity_duration = get(p, "immunity-duration", nothing),
					seed = get(p, "seed", nothing)
				)
				_write_results(res, p["out"])
				return nothing
			end

		#	Schema: MULTI
			if mode == "multi"
				s = ArgParseSettings()
				@add_arg_table! s begin
					"--graphml"
						arg_type = String
						required = true
					"--sas-transformation"
						action = :store_true
					"--infected"
						arg_type = Int
						nargs = '+'
						required = true
					"--inf-r"
						arg_type = Float64
						required = true
					"--rec-t"
						arg_type = Int
						required = true
					"--max-time"
						arg_type = Int
						required = true
					"--p-symp"
						arg_type = Float64
						required = true
					"--b-int"
						arg_type = Float64
						required = true
					"--b-close"
						arg_type = Float64
						required = true
					"--b-cxn-peers"
						arg_type = Float64
						required = true
					"--b-cxn-total"
						arg_type = Float64
						required = true
					"--b-cxn-symp"
						arg_type = Float64
						required = true
					"--b-cls-x-smp"
						arg_type = Float64
						required = true
					"--transmission"
						arg_type = String
						default = "weighted"
					"--immunity-duration"
						arg_type = Int
					"--seed"
						arg_type = Int
					"--out"
						arg_type = String
						default = "-"
				end

				p = parse_args(args[2:end], s)
				if p === nothing; return nothing; end

				graphml_files = String.(strip.(split(p["graphml"], ',')))
				for f in graphml_files
					if !isfile(f)
						error("GraphML file not found: $f")
					end
				end

				res = run_multi(
					graphml_files;
					sas_transformation = get(p, "sas-transformation", false),
					infectedp = Vector{Int}(p["infected"]),
					inf_r = p["inf-r"],
					rec_t = p["rec-t"],
					maxtime = p["max-time"],
					p_symp = p["p-symp"],
					b_int = p["b-int"],
					b_close = p["b-close"],
					b_cxn_peers = p["b-cxn-peers"],
					b_cxn_total = p["b-cxn-total"],
					b_cxn_symp = p["b-cxn-symp"],
					b_cls_x_smp = p["b-cls-x-smp"],
					transmission_method = p["transmission"] == "simple" ? :simple : :weighted,
					immunity_duration = get(p, "immunity-duration", nothing),
					seed = get(p, "seed", nothing)
				)
				_write_results(res, p["out"])
				return nothing
			end

		#	Fallback: unknown mode
			error("Unknown mode: $mode. Use 'single', 'sweep', or 'multi'.")
	end

end # module
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
		#	Extract components
			infection_log = payload["infection_log"]
			metadata = Dict(
				"total_time" => payload["total_time"],
				"final_state_summary" => Dict(
					"n_susceptible" => sum(payload["final_state"][:, 3]),
					"n_infected"    => sum(payload["final_state"][:, 2]),
					"n_recovered"   => sum(payload["final_state"][:, 4])
				)
			)

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

#	Single Wrapper: sim_prep → sirdif pipeline (GraphML in, Dict out)
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

		#	Prepare network (internal)
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

		#	Return Dict (no I/O here)
			return res
	end

#	Sweep Wrapper: vector inputs for coefficients → stacked DataFrame
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

		#	Optional global seed (per-run RNG still varies by iteration order)
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

		#	Triple loop over coefficient grids + iterations
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

		#	Return stacked results + final_state from last run (for metadata summary)
			return Dict(
				"infection_log" => results,
				"total_time"    => maxtime,
				"final_state"   => last_res === nothing ? zeros(Int, 0, 0) : last_res["final_state"]
			)
	end

#	CLI Main Function (modes: single | sweep)
	function cli_main(args::Vector{String} = ARGS)
		#	No args: print short usage
			if isempty(args)
				println("Usage:")
				println("  single run: diffusion_sim.CLI single --graphml <file> ... --out <path>")
				println("  sweep  run: diffusion_sim.CLI sweep  --graphml <file> ... --out <path>")
				return nothing
			end

		#	Mode select
			mode = args[1]

		#	Schema: SINGLE
			if mode == "single"
				s = ArgParseSettings()
				@add_arg_table! s begin
					"--graphml"
						help     = "Path to GraphML file (exported from igraph)."
						arg_type = String
						required = true

					"--sas-transformation"
						help   = "Apply SAS-style preprocessing to the graph."
						action = :store_true

					"--infected"
						help     = "One or more initial infected node IDs."
						arg_type = Int
						nargs    = '+'
						required = true

					"--inf-r";               arg_type = Float64; required = true
					"--rec-t";               arg_type = Int;     required = true
					"--max-time";            arg_type = Int;     required = true
					"--p-symp";              arg_type = Float64; required = true
					"--b-int";               arg_type = Float64; required = true
					"--b-close";             arg_type = Float64; required = true
					"--b-cxn-peers";         arg_type = Float64; required = true
					"--b-cxn-total";         arg_type = Float64; required = true
					"--b-cxn-symp";          arg_type = Float64; required = true
					"--b-cls-x-smp";         arg_type = Float64; required = true

					"--transmission"
						help     = "weighted | simple (default: weighted)"
						arg_type = String
						default  = "weighted"

					"--immunity-duration";   arg_type = Int
					"--seed";                arg_type = Int

					"--out"
						help     = "Output base path ('-' = stdout CSV only)."
						arg_type = String
						default  = "-"
				end

				p = parse_args(args[2:end], s)
				if p === nothing
					return nothing
				end

				#	Convert arguments
					infectedp = Vector{Int}(p["infected"])
					transmission_method = p["transmission"] == "simple" ? :simple : :weighted

				#	Run + write
					res = run_simulation(
						p["graphml"];
						sas_transformation = get(p, "sas-transformation", false),
						infectedp = infectedp,
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
						transmission_method = transmission_method,
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
						help     = "Path to GraphML file (exported from igraph)."
						arg_type = String
						required = true

					"--sas-transformation"
						help   = "Apply SAS-style preprocessing to the graph."
						action = :store_true

					"--infected"
						help     = "One or more initial infected node IDs."
						arg_type = Int
						nargs    = '+'
						required = true

					"--inf-r";               arg_type = Float64; required = true
					"--rec-t";               arg_type = Int;     required = true
					"--max-time";            arg_type = Int;     required = true
					"--p-symp";              arg_type = Float64; required = true
					"--b-int";               arg_type = Float64; required = true
					"--b-close";             arg_type = Float64; required = true

					"--b-cxn-peers"
						help     = "Comma-separated values for peer coefficient (e.g., -0.5,-0.75,-1.0)"
						arg_type = String
						required = true

					"--b-cxn-total"
						help     = "Comma-separated values for global coefficient (e.g., -3.5,-5.0,-6.5)"
						arg_type = String
						required = true

					"--b-cxn-symp"
						help     = "Comma-separated values for symptomatic coefficient (e.g., -1.5,-2.5)"
						arg_type = String
						required = true

					"--b-cls-x-smp";         arg_type = Float64; required = true

					"--num-iter"
						help     = "Number of iterations per (bp,bg,bs) combination."
						arg_type = Int
						default  = 1

					"--transmission"
						help     = "weighted | simple (default: weighted)"
						arg_type = String
						default  = "weighted"

					"--immunity-duration";   arg_type = Int
					"--seed";                arg_type = Int

					"--out"
						help     = "Output base path ('-' = stdout CSV only)."
						arg_type = String
						default  = "-"
				end

				p = parse_args(args[2:end], s)
				if p === nothing
					return nothing
				end

				#	Convert arguments
					infectedp = Vector{Int}(p["infected"])
					transmission_method = p["transmission"] == "simple" ? :simple : :weighted
					bp = parse.(Float64, strip.(split(p["b-cxn-peers"],  ',')))
					bg = parse.(Float64, strip.(split(p["b-cxn-total"],  ',')))
					bs = parse.(Float64, strip.(split(p["b-cxn-symp"],   ',')))

				#	Run sweep + write
					res = run_sweep(
						p["graphml"];
						sas_transformation = get(p, "sas-transformation", false),
						infectedp = infectedp,
						inf_r = p["inf-r"],
						rec_t = p["rec-t"],
						maxtime = p["max-time"],
						p_symp = p["p-symp"],
						b_int = p["b-int"],
						b_close = p["b-close"],
						b_cxn_peers = bp,
						b_cxn_total = bg,
						b_cxn_symp  = bs,
						b_cls_x_smp = p["b-cls-x-smp"],
						num_iter = p["num-iter"],
						transmission_method = transmission_method,
						immunity_duration = get(p, "immunity-duration", nothing),
						seed = get(p, "seed", nothing)
					)
					_write_results(res, p["out"])
					return nothing
			end

		#	Fallback: unknown mode
			error("Unknown mode: $mode. Use 'single' or 'sweep'.")
	end

end # module

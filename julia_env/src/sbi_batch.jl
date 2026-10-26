#Batch Simulation for Simulation-Based Inference
#Jonathan H. Morgan, Ph.D.
#26 October 2026
#Included in module diffusion_sim after richards.jl and summarize.jl. Runs a batch
#of parameter draws on one network, summarizes each run, and writes one JSON line
#per run (schema version 3) with a manifest of the batch settings. This is the
#interface to the Python layer; the single, sweep, and multi CLI modes, which write
#full trajectories as CSV, are unchanged.
#Revised 26 October 2026: sas_transformation option, so that batches can use the
#SAS-preprocessed network (symmetrized, with common-third ties) as the validated
#replication did. Seeds are drawn from the largest weakly connected component by
#default (seed_component), and the manifest records the component structure.

#   Modules
	using .richards_model
	using .run_summary

# ============================================================
#	CONSTANTS
# ============================================================

#	Record Layout (schema version 3)
	const SBI_SCHEMA_VERSION = 3
	const THETA_FIELDS = (:inf_r, :b_int, :b_close, :b_cxn_peers, :b_cxn_total, :b_cxn_symp,
						  :b_cls_x_smp, :p_symp, :rec_scale, :rec_shape, :imm_scale, :imm_shape)
	const SUMMARY_FIELDS = (:K, :r_growth, :t_i, :a, :frac_at_inflection, :fit_rmse,
							:peak_incidence, :peak_prevalence, :peak_infected, :peak_day,
							:infected_at_surge, :max_slope, :final_size, :duration, :extinguished,
							:t_10, :t_half, :t_90, :C0, :fit_ok, :extinct, :viable, :status, :class)
	const TARGET_FIELDS = (:in_target, :peak_ok, :surge_ok, :slope_ok, :prevalence_ok, :final_size_ok)

# ============================================================
#	SEEDS
# ============================================================

#	Helper Function for run_seeds: SplitMix64 mixing step
	function splitmix64(x::UInt64)::UInt64
		"""
		Args:
			x::UInt64: input state
		Returns:
			UInt64: the SplitMix64 output for state x
		Notes:
			Steele, Lea, and Flood (2014). Arithmetic wraps modulo 2^64.
			splitmix64(0) = 0xe220a8397b1dcdaf, the published first output
			for state 0. A fixed integer function is used rather than
			Base.hash, whose output may change between Julia versions.
		"""

		#	Add the golden-ratio increment and mix
			z = x + 0x9e3779b97f4a7c15
			z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
			z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
			return z ⊻ (z >> 31)
	end

#	Run Seeds from the Batch Master Seed
	function run_seeds(master_seed::Integer, draw_id::Integer, replicate::Integer=1)
		"""
		Args:
			master_seed::Integer: batch master seed
			draw_id::Integer: identifier of the parameter draw
			replicate::Integer: replicate number, from 1 (default = 1)
		Returns:
			NamedTuple: (run_seed::Int64, node_seed::Int64), each in [0, 2^53)
		Notes:
			The three inputs are mixed in turn with SplitMix64. run_seed
			drives sirdif through Xoshiro(run_seed); node_seed drives the
			choice of seed nodes through a separate Xoshiro(node_seed). With
			separate streams, a stored run replays exactly as
			sirdif(...; rng=Xoshiro(run_seed)) with its stored seed nodes.
			Both seeds keep the top 53 bits of the mixed value, so they are
			exact in JSON readers that store numbers as doubles (R,
			JavaScript). Replicate 1 is the run an ordinary batch produces.
		"""

		#	Mix the master seed, the draw, and the replicate
			h = splitmix64(master_seed % UInt64)
			h = splitmix64(h ⊻ (draw_id % UInt64))
			h = splitmix64(h ⊻ (replicate % UInt64))

		#	Assembling Result
			return (run_seed = Int64(h >> 11), node_seed = Int64(splitmix64(h) >> 11))
	end
	@doc raw"""
	**Description**
	Derive a run's random number seeds from the batch master seed, the draw, and the replicate.

	**Usage**
	`run_seeds(master_seed, draw_id, replicate=1)`

	**Details**
	The inputs are mixed in turn with SplitMix64, a fixed integer function, so the seeds do not depend on the Julia version, the thread count, or the order in which runs finish. `run_seed` drives the simulation and `node_seed` drives the choice of seed nodes, through separate `Xoshiro` generators. Both are below ``2^{53}`` so that they survive JSON readers that store numbers as doubles.

	**Value**
	`NamedTuple` `(run_seed, node_seed)`.

	**See Also**
	`run_sbi_batch`, `draw_seed_nodes`
	""" run_seeds

#	Helper Function for weak_components: root of a node's set, with path halving
	function find_root!(parent::Vector{Int}, i::Int)::Int
		"""
		Args:
			parent::Vector{Int}: union-find parent array (modified: paths are shortened)
			i::Int: row index
		Returns:
			Int: the root row of i's set
		"""

		#	Walk to the root, halving the path
			while parent[i] != i
				parent[i] = parent[parent[i]]
				i = parent[i]
			end
			return i
	end

#	Weakly Connected Components of a Network
	function weak_components(alst::Matrix{Int})
		"""
		Args:
			alst::Matrix{Int}: adjacency list from sim_prep (ID in column 1, neighbor IDs after, 0 padding)
		Returns:
			NamedTuple: (label::Vector{Int}, sizes::Dict{Int,Int}, largest::Int)
				label[i] is the component label of row i (the root row of its set);
				sizes maps each label to its number of nodes; largest is the label
				of the largest component
		Notes:
			Ties are treated as undirected, so a one-directional tie joins two
			nodes (weak connectivity). Ties to oneself and to IDs not in the
			network are ignored. If two components share the largest size, the
			one containing the earliest row is chosen, so the choice is
			deterministic. Added 26 October 2026.
		"""

		#	Union-find over the rows
			n = size(alst, 1)
			row_of = Dict(alst[i, 1] => i for i in 1:n)
			parent = collect(1:n)
			for i in 1:n, nbr in @view(alst[i, 2:end])
				if nbr > 0 && nbr != alst[i, 1] && haskey(row_of, nbr)
					ri = find_root!(parent, i)
					rj = find_root!(parent, row_of[nbr])
					if ri != rj
						parent[max(ri, rj)] = min(ri, rj)
					end
				end
			end

		#	Labels and sizes
			label = [find_root!(parent, i) for i in 1:n]
			sizes = Dict{Int, Int}()
			for l in label
				sizes[l] = get(sizes, l, 0) + 1
			end

		#	Largest component, earliest row on a tie
			largest = label[1]
			for i in 1:n
				if sizes[label[i]] > sizes[largest]
					largest = label[i]
				end
			end

		#	Assembling Result
			return (label = label, sizes = sizes, largest = largest)
	end
	@doc raw"""
	**Description**
	Weakly connected components of a network in `sim_prep` format, with the largest identified.

	**Usage**
	`weak_components(alst)`

	**Details**
	Ties are treated as undirected. A component found this way is a set of nodes that an infection could in principle connect; a seed outside the largest component is confined to a small group with no further connection.

	**Value**
	`NamedTuple` `(label, sizes, largest)`: a component label for each row, the size of each component, and the label of the largest.

	**See Also**
	`eligible_seed_nodes`
	""" weak_components

#	Nodes Eligible as Seeds
	function eligible_seed_nodes(alst::Matrix{Int}; component::Symbol=:largest)::Vector{Int}
		"""
		Args:
			alst::Matrix{Int}: adjacency list from sim_prep (ID in column 1, neighbor IDs after, 0 padding)
			component::Symbol: :largest (default) or :any
		Returns:
			Vector{Int}: IDs of eligible nodes, in row order
		Notes:
			A node is eligible if it names at least one other node, so that
			as a seed it can transmit; a tie to oneself does not count.
			With component = :largest (26 October 2026), it must also lie in
			the largest weakly connected component. A seed in a small
			component of a few people cannot start an epidemic whatever the
			parameters, so its run says nothing about them; on the synthetic
			village networks about 30% of nodes lie outside the largest
			component. component = :any keeps the earlier rule (any node
			with a tie), for comparison.
		"""

		#	Validation
			if !(component in (:largest, :any))
				throw(ArgumentError("component must be :largest or :any"))
			end

		#	Component of each row, when restricting to the largest
			comps = component == :largest ? weak_components(alst) : nothing

		#	Keep rows with a tie to another node (and in the largest component)
			ids = Int[]
			for i in 1:size(alst, 1)
				ego_id = alst[i, 1]
				has_tie = any(nbr -> nbr > 0 && nbr != ego_id, @view(alst[i, 2:end]))
				in_scope = comps === nothing || comps.label[i] == comps.largest
				if has_tie && in_scope
					push!(ids, ego_id)
				end
			end

		#	Assembling Result
			return ids
	end
	@doc raw"""
	**Description**
	The IDs of nodes that may be drawn as seeds: by default, nodes in the largest weakly connected component that name at least one other node.

	**Usage**
	`eligible_seed_nodes(alst; component=:largest)`

	**Details**
	`component = :largest` (the default since 26 October 2026) excludes seeds confined to small components, where an introduction dies out for structural reasons alone. `component = :any` keeps the earlier rule, any node with a tie to another node.

	**Value**
	`Vector{Int}` of node IDs.

	**See Also**
	`weak_components`, `draw_seed_nodes`, `run_sbi_batch`
	""" eligible_seed_nodes

#	Seed Node Draw
	function draw_seed_nodes(rng::AbstractRNG, eligible::Vector{Int}, n_seed_nodes::Int)::Vector{Int}
		"""
		Args:
			rng::AbstractRNG: generator for the draw (Xoshiro(node_seed) in a batch)
			eligible::Vector{Int}: IDs from eligible_seed_nodes
			n_seed_nodes::Int: number of seed nodes (≥ 1)
		Returns:
			Vector{Int}: distinct seed node IDs
		Notes:
			Drawn without replacement with StatsBase.sample. Throws
			ArgumentError if n_seed_nodes is below 1 or exceeds the number of
			eligible nodes.
		"""

		#	Validation
			if n_seed_nodes < 1 || n_seed_nodes > length(eligible)
				throw(ArgumentError("n_seed_nodes must be between 1 and $(length(eligible)), the number of eligible nodes"))
			end

		#	Draw without replacement
			return sample(rng, eligible, n_seed_nodes; replace=false)
	end
	@doc raw"""
	**Description**
	Draw distinct seed nodes from the eligible nodes.

	**Usage**
	`draw_seed_nodes(rng, eligible, n_seed_nodes)`

	**Value**
	`Vector{Int}` of node IDs.

	**See Also**
	`eligible_seed_nodes`, `run_seeds`
	""" draw_seed_nodes

#	Read a Seed File (testing mode)
	function read_seed_file(seed_file::String)
		"""
		Args:
			seed_file::String: JSON array of objects with draw_id, run_seed,
							   seed_nodes, and optionally replicate (default 1)
		Returns:
			Dict{Tuple{Int,Int}, NamedTuple}: (draw_id, replicate) → (run_seed, seed_nodes)
		Notes:
			Testing mode reads every run's seed and seed nodes from this file,
			so a batch is replicable without a master seed. write_seed_file
			produces such a file from a batch output, which makes any batch
			replayable in testing mode. Throws ArgumentError for a missing
			field or a repeated (draw_id, replicate) pair.
		"""

		#	Parse the file
			entries = JSON3.read(read(seed_file, String))
			table = Dict{Tuple{Int,Int}, NamedTuple{(:run_seed, :seed_nodes), Tuple{Int64, Vector{Int}}}}()

		#	Validate and store each entry
			for (k, e) in enumerate(entries)
				for field in (:draw_id, :run_seed, :seed_nodes)
					if !haskey(e, field)
						throw(ArgumentError("seed file entry $(k) lacks $(field)"))
					end
				end
				key = (Int(e[:draw_id]), haskey(e, :replicate) ? Int(e[:replicate]) : 1)
				if haskey(table, key)
					throw(ArgumentError("seed file repeats draw_id $(key[1]), replicate $(key[2])"))
				end
				table[key] = (run_seed = Int64(e[:run_seed]), seed_nodes = Int.(collect(e[:seed_nodes])))
			end

		#	Assembling Result
			return table
	end
	@doc raw"""
	**Description**
	Read the per-run seeds and seed nodes for a batch in testing mode.

	**Usage**
	`read_seed_file(seed_file)`

	**Details**
	The file is a JSON array; each entry gives `draw_id`, `run_seed`, `seed_nodes`, and optionally `replicate` (default 1).

	**Value**
	`Dict` keyed by `(draw_id, replicate)`.

	**See Also**
	`write_seed_file`, `run_sbi_batch`
	""" read_seed_file

#	Write a Seed File from a Batch Output
	function write_seed_file(out_file::String, seed_file::String)::Int
		"""
		Args:
			out_file::String: JSONL output of run_sbi_batch
			seed_file::String: path of the seed file to write
		Returns:
			Int: number of entries written
		Notes:
			Copies draw_id, replicate, run_seed, and seed_nodes from each
			record, so the batch can be rerun exactly in testing mode.
		"""

		#	Collect the seeds of every record
			entries = []
			for line in eachline(out_file)
				meta = JSON3.read(line)[:meta]
				push!(entries, (draw_id = meta[:draw_id], replicate = meta[:replicate],
								run_seed = meta[:run_seed], seed_nodes = collect(meta[:seed_nodes])))
			end

		#	Write the file
			io = open(seed_file, "w")
			try
				JSON3.write(io, entries)
			finally
				close(io)
			end

		#	Assembling Result
			return length(entries)
	end
	@doc raw"""
	**Description**
	Write a seed file from a batch output, so the batch can be rerun exactly in testing mode.

	**Usage**
	`write_seed_file(out_file, seed_file)`

	**Value**
	`Int` number of entries.

	**See Also**
	`read_seed_file`, `run_sbi_batch`
	""" write_seed_file

# ============================================================
#	PARAMETER DRAWS
# ============================================================

#	Helper Function for read_draws: parse and validate one draw
	function parse_draw(obj, k::Int)
		"""
		Args:
			obj: one JSON object from the parameter file
			k::Int: position of the draw in the file, for error messages
		Returns:
			NamedTuple: (draw_id::Int, theta::NamedTuple) with the 12 parameters as Float64,
						and imm_scale and imm_shape as nothing for permanent immunity
		Notes:
			imm_scale and imm_shape must both be numbers or both be null or
			absent. Throws ArgumentError for a missing field, a non-integer
			draw_id, or a value outside its range. The duration types check
			their own parameters when the run is built.
		"""

		#	Required fields
			for field in (:draw_id, THETA_FIELDS[1:10]...)
				if !haskey(obj, field) || obj[field] === nothing
					throw(ArgumentError("parameter draw $(k) lacks $(field)"))
				end
			end
			if !(obj[:draw_id] isa Integer)
				throw(ArgumentError("parameter draw $(k): draw_id must be an integer"))
			end

		#	Immunity: both numbers, or both null or absent (permanent)
			imm_scale = haskey(obj, :imm_scale) ? obj[:imm_scale] : nothing
			imm_shape = haskey(obj, :imm_shape) ? obj[:imm_shape] : nothing
			if (imm_scale === nothing) != (imm_shape === nothing)
				throw(ArgumentError("parameter draw $(k): imm_scale and imm_shape must both be set or both be null"))
			end

		#	Convert to Float64
			theta = (inf_r = Float64(obj[:inf_r]), b_int = Float64(obj[:b_int]), b_close = Float64(obj[:b_close]),
					 b_cxn_peers = Float64(obj[:b_cxn_peers]), b_cxn_total = Float64(obj[:b_cxn_total]),
					 b_cxn_symp = Float64(obj[:b_cxn_symp]), b_cls_x_smp = Float64(obj[:b_cls_x_smp]),
					 p_symp = Float64(obj[:p_symp]), rec_scale = Float64(obj[:rec_scale]),
					 rec_shape = Float64(obj[:rec_shape]),
					 imm_scale = imm_scale === nothing ? nothing : Float64(imm_scale),
					 imm_shape = imm_shape === nothing ? nothing : Float64(imm_shape))

		#	Probabilities
			if !(0.0 <= theta.inf_r <= 1.0) || !(0.0 <= theta.p_symp <= 1.0)
				throw(ArgumentError("parameter draw $(k): inf_r and p_symp must lie in [0, 1]"))
			end

		#	Assembling Result
			return (draw_id = Int(obj[:draw_id]), theta = theta)
	end

#	Read a Parameter File
	function read_draws(params_file::String)
		"""
		Args:
			params_file::String: JSON array of parameter draws, each with draw_id and the 12 parameters
		Returns:
			Vector{NamedTuple}: parsed draws in file order
		Notes:
			Throws ArgumentError for an empty file or a repeated draw_id.
		"""

		#	Parse every draw
			raw = JSON3.read(read(params_file, String))
			draws = [parse_draw(obj, k) for (k, obj) in enumerate(raw)]

		#	Validation
			if isempty(draws)
				throw(ArgumentError("parameter file contains no draws"))
			end
			ids = [d.draw_id for d in draws]
			if length(unique(ids)) != length(ids)
				throw(ArgumentError("draw_id values must be unique"))
			end

		#	Assembling Result
			return draws
	end
	@doc raw"""
	**Description**
	Read and validate a batch's parameter draws.

	**Usage**
	`read_draws(params_file)`

	**Details**
	The file is a JSON array of objects with `draw_id` and the parameters `inf_r`, `b_int`, `b_close`, `b_cxn_peers`, `b_cxn_total`, `b_cxn_symp`, `b_cls_x_smp`, `p_symp`, `rec_scale`, `rec_shape`, `imm_scale`, and `imm_shape`. The immunity parameters may both be `null` for permanent immunity.

	**Value**
	`Vector` of `(draw_id, theta)` named tuples.

	**See Also**
	`run_sbi_batch`
	""" read_draws

# ============================================================
#	ONE RUN
# ============================================================

#	Simulate and Summarize One Run
	function simulate_and_summarize(alst::Matrix{Int}, vlst::Matrix{Float64}, seed_nodes::Vector{Int},
									theta::NamedTuple, maxtime::Int, region::TargetRegion;
									rng::AbstractRNG,
									transmission_method::Symbol=:weighted,
									extinct_policy::Symbol=:exclude,
									extinction_threshold::Float64=0.02)
		"""
		Args:
			alst::Matrix{Int}, vlst::Matrix{Float64}: network from sim_prep
			seed_nodes::Vector{Int}: IDs of the initially infected nodes
			theta::NamedTuple: the 12 parameters (imm_scale and imm_shape nothing for permanent immunity)
			maxtime::Int: simulation horizon in days
			region::TargetRegion: acceptable outcome region
			rng::AbstractRNG: generator for the simulation (Xoshiro(run_seed) in a batch)
			transmission_method::Symbol: :weighted or :simple (default = :weighted)
			extinct_policy::Symbol: :exclude, :success, or :failure (default = :exclude)
			extinction_threshold::Float64: passed to summarize_run (default = 0.02)
		Returns:
			NamedTuple: the summarize_run result
		Notes:
			Recovery is always Weibull, as the prior draws it. Calls the
			validated sirdif without change. The parameter order follows
			sirdif: b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp,
			b_cls_x_smp.
		"""

		#	Duration specifications
			recovery = WeibullDuration(theta.rec_scale, theta.rec_shape)
			immunity = theta.imm_scale === nothing ? nothing : WeibullDuration(theta.imm_scale, theta.imm_shape)

		#	Simulate
			res = sirdif(alst, vlst, seed_nodes, theta.inf_r, recovery, maxtime, theta.p_symp,
						 theta.b_int, theta.b_close, theta.b_cxn_peers, theta.b_cxn_total,
						 theta.b_cxn_symp, theta.b_cls_x_smp;
						 transmission_method=transmission_method, immunity=immunity, rng=rng)

		#	Summarize
			return summarize_run(res["infection_log"], size(alst, 1), maxtime, region;
								 rec_scale=theta.rec_scale, rec_shape=theta.rec_shape,
								 extinct_policy=extinct_policy, extinction_threshold=extinction_threshold)
	end
	@doc raw"""
	**Description**
	Run one simulation from a parameter draw and summarize it: `sirdif`, then `summarize_run`.

	**Usage**
	```julia
	simulate_and_summarize(alst, vlst, seed_nodes, theta, maxtime, region;
						   rng, transmission_method=:weighted,
						   extinct_policy=:exclude, extinction_threshold=0.02)
	```

	**Details**
	To replay a stored run, pass its `seed_nodes` and `rng=Xoshiro(run_seed)` with the draw's parameters and the batch settings from the manifest.

	**Value**
	The `NamedTuple` returned by `summarize_run`.

	**See Also**
	`run_sbi_batch`, `summarize_run`, `sirdif`
	""" simulate_and_summarize

#	Helper Function for build_record: one JSON value
	function json_value(x)
		"""
		Args:
			x: a value from a summary or parameter tuple
		Returns:
			The value to write: nothing for a non-finite float (JSON null),
			a String for a Symbol, and the value itself otherwise
		Notes:
			JSON has no NaN or Inf, and JSON3 refuses to write them.
		"""

		#	Convert by type
			if x isa AbstractFloat
				return isfinite(x) ? x : nothing
			elseif x isa Symbol
				return String(x)
			end
			return x
	end

#	Helper Function for build_record: one block of named fields
	function json_block(nt::NamedTuple, fields::Tuple)
		"""
		Args:
			nt::NamedTuple: source tuple
			fields::Tuple: field names to keep, in output order
		Returns:
			NamedTuple: the selected fields, each passed through json_value
		"""

		#	Select and convert
			return NamedTuple{fields}(Tuple(json_value(getfield(nt, f)) for f in fields))
	end

#	Helper Function for run_sbi_batch: assemble one output record
	function build_record(theta::NamedTuple, summary::NamedTuple, meta::NamedTuple)
		"""
		Args:
			theta::NamedTuple: the parsed parameters of the draw
			summary::NamedTuple: result of simulate_and_summarize
			meta::NamedTuple: village, draw_id, replicate, run_seed, seed_nodes, n_nodes, schema_version
		Returns:
			NamedTuple: (theta, summary, target, meta), ready for JSON3.write
		"""

		#	Assembling Result
			return (theta = json_block(theta, THETA_FIELDS),
					summary = json_block(summary, SUMMARY_FIELDS),
					target = json_block(summary, TARGET_FIELDS),
					meta = meta)
	end

#	Helper Function for run_sbi_batch: simulate one task and return its JSON line
	function run_task(draw::NamedTuple, replicate::Int, batch::NamedTuple)
		"""
		Args:
			draw::NamedTuple: (draw_id, theta) from read_draws
			replicate::Int: replicate number
			batch::NamedTuple: shared, read-only batch settings and network
		Returns:
			Tuple{String, Symbol}: the JSON line and the run class
		Notes:
			Called concurrently from several threads. It reads only the
			shared settings and builds its own generators, so runs do not
			interact.
		"""

		#	Seeds for this run
			if batch.seed_mode == :random
				s = run_seeds(batch.master_seed, draw.draw_id, replicate)
				run_seed = s.run_seed
				seed_nodes = draw_seed_nodes(Xoshiro(s.node_seed), batch.eligible, batch.n_seed_nodes)
			else
				entry = batch.seed_table[(draw.draw_id, replicate)]
				run_seed = entry.run_seed
				seed_nodes = entry.seed_nodes
			end

		#	Simulate and summarize
			summary = simulate_and_summarize(batch.alst, batch.vlst, seed_nodes, draw.theta, batch.maxtime, batch.region;
											 rng=Xoshiro(run_seed), transmission_method=batch.transmission_method,
											 extinct_policy=batch.extinct_policy,
											 extinction_threshold=batch.extinction_threshold)

		#	Build the record
			meta = (village = batch.village, draw_id = draw.draw_id, replicate = replicate, run_seed = run_seed,
					seed_nodes = seed_nodes, n_nodes = batch.n_nodes, schema_version = SBI_SCHEMA_VERSION)
			record = build_record(draw.theta, summary, meta)

		#	Assembling Result
			return (JSON3.write(record), summary.class)
	end

#	Helper Function for run_sbi_batch: manifest path beside the output
	function manifest_path(out_file::String)::String
		"""
		Args:
			out_file::String: path of the JSONL output
		Returns:
			String: results.jsonl → results.manifest.json; any other name gains .manifest.json
		"""

		#	Replace the extension
			return endswith(out_file, ".jsonl") ? out_file[1:end-6] * ".manifest.json" : out_file * ".manifest.json"
	end

# ============================================================
#	BATCH
# ============================================================

#	Run a Batch of Parameter Draws
	function run_sbi_batch(graphml_file::String, params_file::String, out_file::String;
						   village::String,
						   sas_transformation::Bool=false,
						   seed_component::Symbol=:largest,
						   region::TargetRegion=TargetRegion(),
						   maxtime::Int=200,
						   transmission_method::Symbol=:weighted,
						   seed_mode::Symbol=:random,
						   master_seed::Union{Integer, Nothing}=nothing,
						   seed_file::Union{String, Nothing}=nothing,
						   n_seed_nodes::Int=1,
						   replicates::Int=1,
						   extinct_policy::Symbol=:exclude,
						   extinction_threshold::Float64=0.02,
						   threaded::Bool=true)
		"""
		Args:
			graphml_file::String: network file
			params_file::String: JSON array of parameter draws (read_draws)
			out_file::String: JSONL output, one line per run
			village::String: village identifier for the records and manifest
			sas_transformation::Bool: load the network with SAS preprocessing (default = false)
			seed_component::Symbol: :largest (default) draws seeds from the largest weakly
									connected component; :any from any node with a tie
			region::TargetRegion: acceptable outcome region (default = unconstrained)
			maxtime::Int: simulation horizon in days (default = 200)
			transmission_method::Symbol: :weighted or :simple (default = :weighted)
			seed_mode::Symbol: :random (master seed) or :file (seed file) (default = :random)
			master_seed: batch master seed, required in random mode
			seed_file: seed file, required in file mode
			n_seed_nodes::Int: seed nodes per run in random mode (default = 1)
			replicates::Int: runs per draw (default = 1; SBI batches use 1)
			extinct_policy::Symbol: :exclude, :success, or :failure (default = :exclude)
			extinction_threshold::Float64: growth below which a run is extinct (default = 0.02)
			threaded::Bool: spread runs across Julia's threads (default = true)
		Returns:
			NamedTuple: (out_file, manifest_file, n_runs, class_counts::Dict{Symbol,Int})
		Notes:
			Runs are independent, so they are spread across threads; a
			single run is never split, because its dynamics depend on order.
			Each thread writes into its own slot, and the file is written in
			draw order, then replicate order, after every run finishes. The
			output therefore does not depend on the thread count. BLAS is
			held to one thread during the batch so that the curve fits do
			not oversubscribe the CPU; the previous setting is restored.
			sas_transformation passes through to sim_prep. The SAS form
			symmetrizes the network and adds common-third ties with
			square-root weights, and it keeps only nodes that appear in an
			edge, so isolated nodes drop out and N, C0, and the eligible
			seeds follow the transformed network. The manifest records the
			choice, and a replay must load the network the same way.
		"""

		#	Validation of the settings
			if !(seed_mode in (:random, :file))
				throw(ArgumentError("seed_mode must be :random or :file"))
			end
			if seed_mode == :random && master_seed === nothing
				throw(ArgumentError("random seed mode requires master_seed"))
			end
			if seed_mode == :file && seed_file === nothing
				throw(ArgumentError("file seed mode requires seed_file"))
			end
			if !(extinct_policy in (:exclude, :success, :failure))
				throw(ArgumentError("extinct_policy must be :exclude, :success, or :failure"))
			end
			if replicates < 1 || maxtime < 1
				throw(ArgumentError("replicates and maxtime must be at least 1"))
			end

		#	Network, eligible seeds, and draws
			network = sim_prep(graphml_file; sas_transformation=sas_transformation)
			alst = network.alst
			vlst = network.vlst
			n_nodes = size(alst, 1)
			eligible = eligible_seed_nodes(alst; component=seed_component)
			comps = weak_components(alst)
			if seed_mode == :random && !(1 <= n_seed_nodes <= length(eligible))
				throw(ArgumentError("n_seed_nodes must be between 1 and $(length(eligible))"))
			end
			draws = read_draws(params_file)

		#	Seed table for file mode, checked against the batch
			seed_table = nothing
			if seed_mode == :file
				seed_table = read_seed_file(seed_file)
				valid_ids = Set(alst[:, 1])
				for d in draws, rep in 1:replicates
					entry = get(seed_table, (d.draw_id, rep), nothing)
					if entry === nothing
						throw(ArgumentError("seed file lacks draw_id $(d.draw_id), replicate $(rep)"))
					end
					if isempty(entry.seed_nodes) || !all(in(valid_ids), entry.seed_nodes)
						throw(ArgumentError("seed file: invalid seed nodes for draw_id $(d.draw_id)"))
					end
				end
			end

		#	Shared, read-only batch settings
			batch = (alst = alst, vlst = vlst, n_nodes = n_nodes, eligible = eligible, region = region,
					 maxtime = maxtime, transmission_method = transmission_method, seed_mode = seed_mode,
					 master_seed = master_seed, seed_table = seed_table, n_seed_nodes = n_seed_nodes,
					 extinct_policy = extinct_policy, extinction_threshold = extinction_threshold,
					 village = village)

		#	Tasks in output order: draws, then replicates
			tasks = [(k, rep) for k in eachindex(draws) for rep in 1:replicates]
			lines = Vector{String}(undef, length(tasks))
			classes = Vector{Symbol}(undef, length(tasks))

		#	Run every task, holding BLAS to one thread
			blas_threads = BLAS.get_num_threads()
			BLAS.set_num_threads(1)
			try
				if threaded && Threads.nthreads() > 1
					Threads.@threads for j in eachindex(tasks)
						lines[j], classes[j] = run_task(draws[tasks[j][1]], tasks[j][2], batch)
					end
				else
					for j in eachindex(tasks)
						lines[j], classes[j] = run_task(draws[tasks[j][1]], tasks[j][2], batch)
					end
				end
			finally
				BLAS.set_num_threads(blas_threads)
			end

		#	Write the records in task order
			io = open(out_file, "w")
			try
				for line in lines
					write(io, line, '\n')
				end
			finally
				close(io)
			end

		#	Write the manifest
			manifest = (schema_version = SBI_SCHEMA_VERSION, village = village, graphml = graphml_file,
						sas_transformation = sas_transformation, n_nodes = n_nodes,
						n_components = length(comps.sizes), largest_component = comps.sizes[comps.largest],
						seed_component = String(seed_component), n_eligible_seeds = length(eligible),
						maxtime = maxtime, transmission = String(transmission_method),
						seed_mode = String(seed_mode), master_seed = master_seed, seed_file = seed_file,
						n_seed_nodes = n_seed_nodes, replicates = replicates, n_draws = length(draws),
						n_runs = length(tasks), extinct_policy = String(extinct_policy),
						extinction_threshold = extinction_threshold,
						target = (I_peak_max = region.I_peak_max, I_surge_max = region.I_surge_max,
								  surge_time = region.surge_time, S_max = region.S_max,
								  P_max_thresh = region.P_max_thresh, K_max = region.K_max),
						julia_version = string(VERSION))
			manifest_file = manifest_path(out_file)
			io = open(manifest_file, "w")
			try
				JSON3.write(io, manifest)
			finally
				close(io)
			end

		#	Count the run classes
			class_counts = Dict{Symbol, Int}()
			for c in classes
				class_counts[c] = get(class_counts, c, 0) + 1
			end

		#	Assembling Result
			return (out_file = out_file, manifest_file = manifest_file, n_runs = length(tasks),
					class_counts = class_counts)
	end
	@doc raw"""
	**Description**
	Simulate, summarize, and evaluate a batch of parameter draws on one network, writing one JSON line per run (schema version 3) and a manifest of the batch settings. This is the output the Python layer reads.

	**Usage**
	```julia
	run_sbi_batch(graphml_file, params_file, out_file; village, sas_transformation=false,
				  seed_component=:largest,
				  region=TargetRegion(), maxtime=200, transmission_method=:weighted,
				  seed_mode=:random, master_seed=nothing, seed_file=nothing,
				  n_seed_nodes=1, replicates=1, extinct_policy=:exclude,
				  extinction_threshold=0.02, threaded=true)
	```

	**Arguments**
	- `graphml_file`, `params_file`, `out_file`: Network, parameter draws (JSON array), and JSONL output.
	- `village::String`: Identifier written to every record and the manifest.
	- `sas_transformation::Bool`: Load the network with SAS preprocessing, as the validated replication did (default `false`). Isolated nodes drop out of the transformed network.
	- `seed_component::Symbol`: `:largest` (default) draws seeds from the largest weakly connected component; `:any` from any node with a tie. The manifest records the number of components, the size of the largest, and the number of eligible seeds.
	- `region::TargetRegion`: The acceptable outcome region.
	- `seed_mode`: `:random` derives each run's seeds from `master_seed`; `:file` reads them from `seed_file` (testing mode).
	- `n_seed_nodes`: Seed nodes per run, drawn from nodes with at least one tie.
	- `replicates`: Runs per draw. SBI batches use 1; more are for troubleshooting.
	- `extinct_policy`, `extinction_threshold`: Treatment of extinct runs.
	- `threaded`: Spread runs across Julia's threads (start Julia with `-t N`).

	**Details**
	Each record has four blocks: `theta` (the 12 parameters as given), `summary` (Richards fit, shape, timing, trajectory summaries, and run class), `target` (`in_target` and each constraint), and `meta` (village, `draw_id`, `replicate`, `run_seed`, `seed_nodes`, `n_nodes`, `schema_version`). Missing numbers and excluded results are written as `null`. The manifest, `results.manifest.json` beside `results.jsonl`, holds the settings shared by the batch.

	Runs are spread across threads, but each run's randomness comes from its own generators, derived from the master seed, `draw_id`, and replicate. The file is written in draw order after every run finishes, so the output is identical for any thread count. Any run replays exactly with `simulate_and_summarize` and `rng=Xoshiro(run_seed)`.

	**Value**
	`NamedTuple` `(out_file, manifest_file, n_runs, class_counts)`.

	**See Also**
	`simulate_and_summarize`, `run_seeds`, `write_seed_file`, `read_draws`
	""" run_sbi_batch

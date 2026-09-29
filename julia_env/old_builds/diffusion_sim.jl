module diffusion_sim
#   Packages
    using CSV
    using DataFrames
	using EzXML
	using Random
    using Statistics
	using StatsBase

#   SUPPORT FUNCTIONS

#	Helper Function for sim_prep: load and parse GraphML file
	function load_graphml(filepath::String)
		"""
		Args:
			filepath::String: path to .graphml file
		Returns:
			NamedTuple: (nodes, edges, weights, node_labels) extracted from GraphML
		Notes:
			Parses GraphML directly to extract network structure.
			Handles both weighted and unweighted networks.
		"""
		
		#	Read and parse XML
			doc = readxml(filepath)
			root = doc.root
			
		#	Find graph element
			graph = nothing
			for child in eachelement(root)
				if child.name == "graph"
					graph = child
					break
				end
			end
			
			if graph === nothing
				error("No graph element found in GraphML file")
			end
			
		#	Determine if directed
			is_directed = haskey(graph, "edgedefault") && graph["edgedefault"] == "undirected" ? false : true
			
		#	Extract nodes by direct traversal
			nodes_xml = []
			node_ids = String[]
			node_mapping = Dict{String, Int}()
			
			for elem in eachelement(graph)
				if elem.name == "node"
					push!(nodes_xml, elem)
					node_id = elem["id"]
					push!(node_ids, node_id)
					node_mapping[node_id] = length(node_ids)
				end
			end
			
		#	Extract edges by direct traversal
			edges_xml = []
			el1 = Int[]
			el2 = Int[]
			elwgt = Float64[]
			
			for elem in eachelement(graph)
				if elem.name == "edge"
					source = elem["source"]
					target = elem["target"]
					
					#	Map string IDs to integer indices
						if haskey(node_mapping, source) && haskey(node_mapping, target)
							src_idx = node_mapping[source]
							tgt_idx = node_mapping[target]
							
							push!(el1, src_idx)
							push!(el2, tgt_idx)
							
							#	Look for weight in data elements
								weight = 1.0
								for data_elem in eachelement(elem)
									if data_elem.name == "data" && haskey(data_elem, "key")
										# Check for various weight key names
										if data_elem["key"] in ["e_weight", "weight", "value", "d1"]
											weight_str = nodecontent(data_elem)
											if !isempty(weight_str)
												weight = parse(Float64, weight_str)
											end
											break
										end
									end
								end
								push!(elwgt, weight)
							
							#	Add reverse edge if undirected
								if !is_directed
									push!(el1, tgt_idx)
									push!(el2, src_idx)
									push!(elwgt, weight)
								end
						else
							println("Warning: Edge references non-existent node(s): $source -> $target")
						end
				end
			end
			
		#	Clean up
			EzXML.finalize(doc) 
			
		#	Check if we found data
			if isempty(node_ids)
				error("No nodes found in GraphML file")
			end
			if isempty(el1)
				error("No edges found in GraphML file")
			end
			println("Loaded network: $(length(node_ids)) nodes, $(length(el1)) directed edges")
			
		#	Return extracted data
			return (
				nodes = collect(1:length(node_ids)),
				edges = (el1, el2),
				weights = elwgt,
				node_labels = node_ids
			)
	end

#	Helper Function for sim_prep: convert valued edgelist to adjacency matrix
	function el2adjval(valid::Vector{Int}, el1::Vector{Int}, el2::Vector{Int}, elwgt::Vector{Float64})
		"""
		Args:
			valid::Vector{Int}: unique node identifiers
			el1::Vector{Int}: source nodes of edges
			el2::Vector{Int}: target nodes of edges
			elwgt::Vector{Float64}: edge weights
		Returns:
			Tuple{Matrix{Float64}, Vector{Int}}: weighted adjacency matrix and node ordering
		Notes:
			Converts edge list representation to matrix form.
			Matrix rows/columns ordered by valid node set.
		"""
		
		#	Initialize adjacency matrix
			nodeset = unique(valid)
			n = length(nodeset)
			adjmat = zeros(Float64, n, n)
			
		#	Create node index mapping
			node_to_idx = Dict(node => idx for (idx, node) in enumerate(nodeset))
		
		#	Populate adjacency matrix
			for i in 1:length(el1)
				iv = el1[i]
				jv = el2[i]
				wv = elwgt[i]
				
				#	Map nodes to matrix indices
					if haskey(node_to_idx, iv) && haskey(node_to_idx, jv)
						iloc = node_to_idx[iv]
						jloc = node_to_idx[jv]
						adjmat[iloc, jloc] = wv
					end
			end
		
		#	Return matrix with node ordering
			return (adjmat, nodeset)
	end

#	Helper Function for preprocess_network_sas_style: convert adjacency to matrix
	function adjacency_list_to_matrix(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: adjacency list (node ID in column 1, neighbors in 2:end)
			vlst::Matrix{Float64}: value list (matching structure)
		Returns:
			Tuple{Matrix{Float64}, Vector{Int}}: (adjacency_matrix, node_ids)
		Notes:
			Converts adjacency list format to square adjacency matrix.
		"""
		
		#	Get node IDs from first column
			node_ids = alst[:, 1]
			n = length(node_ids)
			
		#	Initialize adjacency matrix
			adj_matrix = zeros(Float64, n, n)
			
		#	Fill matrix from adjacency list
			for i in 1:n
				for j in 2:size(alst, 2)
					neighbor = alst[i, j]
					if neighbor > 0
						#	Find neighbor's row index
							neighbor_idx = findfirst(x -> x == neighbor, node_ids)
							if neighbor_idx !== nothing
								adj_matrix[i, neighbor_idx] = vlst[i, j]
							end
					end
				end
			end
			
		#	Return matrix and node ordering
			return (adj_matrix, node_ids)
	end

#	Helper Function for preprocess_network_sas_style: convert matrix back to adjacency list
	function matrix_to_adjacency_list(adj_matrix::Matrix{Float64}, node_ids::Vector{Int})
		"""
		Args:
			adj_matrix::Matrix{Float64}: weighted adjacency matrix
			node_ids::Vector{Int}: node IDs corresponding to matrix rows/columns
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (alst, vlst)
		Notes:
			Converts matrix format back to adjacency list format.
		"""
		
		#	Calculate maximum degree
			n = size(adj_matrix, 1)
			max_degree = maximum(sum(adj_matrix .> 0, dims=2))
			
		#	Initialize adjacency and value lists
			alst = zeros(Int, n, Int(max_degree) + 1)
			vlst = zeros(Float64, n, Int(max_degree) + 1)
			
		#	Fill node IDs in first column
			alst[:, 1] = node_ids
			vlst[:, 1] = Float64.(node_ids)
			
		#	Fill adjacency and value lists
			for i in 1:n
				neighbors = findall(x -> x > 0, adj_matrix[i, :])
				if !isempty(neighbors)
					for (col, j) in enumerate(neighbors)
						if col + 1 <= size(alst, 2)
							alst[i, col + 1] = node_ids[j]
							vlst[i, col + 1] = adj_matrix[i, j]
						end
					end
				end
			end
			
		#	Return lists
			return (alst, vlst)
	end

#	Preprocess Network to Match SAS IML Code
	function preprocess_network_sas_style(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: raw adjacency list from sim_prep
			vlst::Matrix{Float64}: raw value list from sim_prep  
		Returns:
			Tuple{Matrix{Int}, Matrix{Float64}}: (processed_alst, processed_vlst)
		Notes:
			Replicates SAS preprocessing:
			1. Symmetrizes the network
			2. Adds common third-party connections
			3. Takes square root of common thirds
			4. Returns in adjacency list format
		"""
		
		#	Convert to matrix format
			adj_matrix, node_ids = adjacency_list_to_matrix(alst, vlst)
			
		#	Symmetrize: symmat = (amat + amat`)
			symmat = adj_matrix + transpose(adj_matrix)
			
		#	Calculate common thirds: com3rds = ((symmat>0)*(symmat>0))#(symmat>0)
		#	This is the element-wise product of the binary adjacency matrix with itself
		#	In SAS: A#B is element-wise multiplication, * is matrix multiplication
			binary_mat = Float64.(symmat .> 0)
			com3rds = (binary_mat * binary_mat) .* binary_mat
			
		#	Take square root of common thirds: com3rds = com3rds##(0.5)
		#	In SAS: ## is element-wise power
			com3rds = sqrt.(com3rds)
			
		#	Add common thirds to symmetrized matrix
			symmat = symmat + com3rds
			
		#	Convert back to adjacency list format
			processed_alst, processed_vlst = matrix_to_adjacency_list(symmat, node_ids)
			
		#	Return processed networks
			return (processed_alst, processed_vlst)
	end
	@doc raw"""
	**Description**
	Preprocesses a network to match the exact transformations performed in the SAS IML code.
	This ensures Julia simulations use the same network structure as SAS simulations.

	**Usage**
	`preprocess_network_sas_style(alst, vlst)`

	**Arguments**
	- `alst::Matrix{Int}`: Raw adjacency list from `sim_prep`
	- `vlst::Matrix{Float64}`: Raw value list from `sim_prep`

	**Details**
	The function performs the following transformations to match SAS:
	
	1. **Symmetrization**: Creates an undirected network by adding the adjacency matrix 
	   to its transpose: `symmat = A + A'`
	
	2. **Common third-party connections**: Identifies paths of length 2 between nodes.
	   For nodes i and j, this counts how many nodes k exist where both i→k and k→j.
	   Computed as: `com3rds = (A_binary * A_binary) ⊙ A_binary`
	
	3. **Square root weighting**: Takes element-wise square root of common thirds
	   to moderate their influence: `com3rds = sqrt(com3rds)`
	
	4. **Final network**: Adds common thirds to symmetrized network:
	   `final = symmat + com3rds`
	
	This preprocessing strengthens connections between nodes that share common neighbors,
	which can significantly affect epidemic dynamics by creating additional transmission
	pathways and stronger clustering.

	**Value**
	Returns a tuple `(processed_alst, processed_vlst)`:
	- `processed_alst::Matrix{Int}`: Preprocessed adjacency list
	- `processed_vlst::Matrix{Float64}`: Preprocessed edge weights

	**Examples**
	```julia
	# Load raw network
	network_data = sim_prep("network.graphml")
	
	# Apply SAS-style preprocessing
	alst_processed, vlst_processed = preprocess_network_sas_style(
	    network_data.alst, 
	    network_data.vlst
	)
	
	# Use in simulation (now matches SAS)
	results = sirdif(alst_processed, vlst_processed, [1], 
	                0.02, 14, 200, 0.75,
	                -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
	
	# Compare network properties
	raw_edges = sum(network_data.vlst .> 0)
	processed_edges = sum(vlst_processed .> 0)
	println("Raw edges: $raw_edges, Processed edges: $processed_edges")
	```

	**See Also**
	`sim_prep`, `sirdif`, `replicate_sas_simulation`
	""" preprocess_network_sas_style

#	SIMULATION FUNCTIONS

#	Network Data Preparation for SIR Simulation from GraphML
	function sim_prep(graphml_file::String; sas_transformation::Bool = false)
		"""
		Args:
			graphml_file::String: path to .graphml file
			sas_transformation::Bool: apply SAS-style preprocessing (default false)
		Returns:
			NamedTuple: (alst=adjacency_list, vlst=value_list)
		Notes:
			Loads GraphML file and converts to adjacency list format for sir_diffusion.
			If sas_transformation=true, applies symmetrization and common thirds.
		"""
		
		#	Load network from GraphML
			network_data = load_graphml(graphml_file)
			nl = network_data.nodes
			el1, el2 = network_data.edges
			elwgt = network_data.weights
			n = length(nl)
		
		#	Apply SAS transformation if requested
			if sas_transformation
				#	CRITICAL: load_graphml already doubles edges for undirected graphs
				#	We need to deduplicate to match what SAS pajread produces
				#	Keep only unique undirected edges
					edge_set = Set{Tuple{Int,Int}}()
					el1_dedup = Int[]
					el2_dedup = Int[]
					
					for (src, tgt) in zip(el1, el2)
						edge_key = src <= tgt ? (src, tgt) : (tgt, src)
						if !(edge_key in edge_set)
							push!(edge_set, edge_key)
							push!(el1_dedup, src)
							push!(el2_dedup, tgt)
						end
					end
				
				#	Build adjacency matrix exactly as SAS adj() does
					unique_nodes = sort(unique(vcat(el1_dedup, el2_dedup)))
					n = length(unique_nodes)
					node_to_idx = Dict(node => idx for (idx, node) in enumerate(unique_nodes))
					
				#	Create adjacency matrix (counting edges)
					amat = zeros(Float64, n, n)
					for (sender, receiver) in zip(el1_dedup, el2_dedup)
						i = node_to_idx[sender]
						j = node_to_idx[receiver]
						amat[i, j] += 1.0  # SAS adj() increments for each edge
					end
				
				#	Step 2a: Symmetrize (amat + amat')
					symmat = amat + transpose(amat)
					
				#	Step 2b: Calculate common thirds
				#	com3rds = ((symmat>0)*(symmat>0))#(symmat>0)
					binary = Float64.(symmat .> 0)
					com3rds = (binary * binary) .* binary
					
				#	Step 2c: Square root of common thirds
				#	com3rds = com3rds##0.5
					com3rds = sqrt.(com3rds)
					
				#	Step 2d: Add to symmetrized matrix
				#	symmat = symmat + com3rds
					symmat = symmat + com3rds
					
				#	Convert to adjacency list format
					max_degree = Int(maximum(sum(symmat .> 0, dims=2)))
					alst = zeros(Int, n, max_degree + 1)
					vlst = zeros(Float64, n, max_degree + 1)
					
				#	Fill ID columns (using the sorted unique nodes)
					alst[:, 1] = unique_nodes
					vlst[:, 1] = Float64.(unique_nodes)
					
				#	Fill neighbor lists
					for i in 1:n
						neighbors = findall(x -> x > 0, symmat[i, :])
						for (col_idx, j) in enumerate(neighbors)
							if col_idx + 1 <= size(alst, 2)
								alst[i, col_idx + 1] = unique_nodes[j]
								vlst[i, col_idx + 1] = symmat[i, j]
							end
						end
					end
					
					println("Applied SAS transformation: symmetrization + common thirds")
			else
				#	No SAS transformation - use raw network
				#	Create binary adjacency matrix for structure
					adj_mat_binary, nodeset = el2adjval(nl, el1, el2, ones(Float64, length(elwgt)))
					
				#	Create weighted adjacency matrix
					adj_mat_weighted, _ = el2adjval(nl, el1, el2, elwgt)
					
				#	Calculate maximum degree
					max_degree = Int(maximum(sum(adj_mat_binary .> 0, dims=2)))
					
				#	Initialize adjacency list and value list
					alst = zeros(Int, n, max_degree + 1)  # +1 for node ID column
					vlst = zeros(Float64, n, max_degree + 1)
					
				#	Populate node IDs in first column
					for i in 1:n
						alst[i, 1] = nodeset[i]
						vlst[i, 1] = Float64(nodeset[i])
					end
				
				#	Fill adjacency and value lists
					for i in 1:n
						#	Find neighbors
							neighbors = findall(x -> x > 0, adj_mat_binary[i, :])
							
							if !isempty(neighbors)
								#	Map indices to node IDs and extract weights
									neighbor_ids = [nodeset[j] for j in neighbors]
									weights = [adj_mat_weighted[i, j] for j in neighbors]
									
								#	Fill adjacency list
									for (col, nid) in enumerate(neighbor_ids)
										if col + 1 <= size(alst, 2)
											alst[i, col + 1] = nid
										end
									end
									
								#	Fill value list
									for (col, wgt) in enumerate(weights)
										if col + 1 <= size(vlst, 2)
											vlst[i, col + 1] = wgt
										end
									end
							end
					end
			end
		
		#	Return formatted lists for simulation
			return (alst = alst, vlst = vlst)
	end
	@doc raw"""
	**Description**
	Prepares network data from a GraphML file for use in SIR diffusion simulations.
	Directly parses GraphML format and converts to adjacency list representation with
	optional SAS-style preprocessing for comparison with SAS PROC IML simulations.

	**Usage**
	`sim_prep(graphml_file; sas_transformation=false)`

	**Arguments**
	- `graphml_file::String`: Path to .graphml file containing network data
	- `sas_transformation::Bool`: Apply SAS-style preprocessing (default `false`)

	**Details**
	The function performs the following operations:
	1. Parses GraphML file using EzXML to extract nodes, edges, and weights
	2. Handles both directed and undirected graphs automatically
	3. If `sas_transformation=true`, applies these transformations:
	- **Symmetrization**: Creates undirected network via `A + A'`
	- **Common thirds**: Adds edges between nodes sharing neighbors via `(A*A) ⊙ A`
	- **Square root weighting**: Applies `sqrt()` to common third connections
	- **Combination**: Final network = symmetrized + sqrt(common thirds)
	4. Converts to adjacency list format where:
	- Each row represents a node
	- First column contains the node ID
	- Subsequent columns contain neighbor IDs (0 for empty slots)
	5. Creates parallel value list containing edge weights

	The SAS transformation significantly changes network topology by adding edges
	between nodes that share common neighbors, which can increase epidemic spread
	by creating additional transmission pathways.

	GraphML weight attributes are detected automatically. The function looks for
	edge attributes named "weight", "value", "e_weight", or "d1". If no weights 
	are found, all edges are assigned weight 1.0.

	**Value**
	Returns a NamedTuple with two fields:
	- `alst`: Matrix{Int} where row i contains node i's ID followed by its neighbor IDs
	- `vlst`: Matrix{Float64} with same structure containing edge weights

	**Examples**
	```julia
	# Load network without SAS preprocessing (raw network)
	network_data = sim_prep("network.graphml")

	# Load network with SAS preprocessing (for comparison with SAS)
	network_sas = sim_prep("network.graphml", sas_transformation=true)

	# Compare network properties
	raw_edges = sum(network_data.vlst .> 0)
	sas_edges = sum(network_sas.vlst .> 0)
	println("Raw edges: $raw_edges, SAS-processed edges: $sas_edges")

	# Use in simulations
	# For Julia-only studies:
	results_julia = sirdif(network_data.alst, network_data.vlst, [1], 
						0.02, 14, 200, 0.75, -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)

	# For SAS comparison studies:
	results_sas_style = sirdif(network_sas.alst, network_sas.vlst, [1], 
							0.02, 14, 200, 0.75, -0.1, 1.0, -0.5, -3.5, -1.5, -0.1)
	```

	**See Also**
	`sirdif`, `replicate_sas_simulation`, `sas_simulation_comparer`
	""" sim_prep

#	Helper Function for sirdif: convert matrix adjacency to vector format
	function matrix_to_vector_adjacency(alst::Matrix{Int}, vlst::Matrix{Float64})
		"""
		Args:
			alst::Matrix{Int}: adjacency matrix from sim_prep
			vlst::Matrix{Float64}: value matrix from sim_prep
		Returns:
			Tuple{Vector{Vector{Int}}, Vector{Vector{Float64}}}: (alst_vec, vlst_vec)
		Notes:
			Converts matrix format to vector format for efficiency.
			Pre-allocates output vectors.
		"""
		
		#	Get dimensions
			n = size(alst, 1)
			
		#	Pre-allocate vector versions
			alst_vec = Vector{Vector{Int}}(undef, n)
			vlst_vec = Vector{Vector{Float64}}(undef, n)
			
		#	Convert each row
			for i in 1:n
				#	Count non-zero neighbors first
					n_neighbors = count(x -> x > 0, @view alst[i, 2:end])
					
				#	Pre-allocate neighbor and weight vectors
					neighbors = Vector{Int}(undef, n_neighbors)
					weights = Vector{Float64}(undef, n_neighbors)
					
				#	Fill pre-allocated vectors
					idx = 1
					for j in 2:size(alst, 2)
						if alst[i, j] > 0
							neighbors[idx] = alst[i, j]
							weights[idx] = vlst[i, j]
							idx += 1
						end
					end
					
				#	Assign to output
					alst_vec[i] = neighbors
					vlst_vec[i] = weights
			end
			
		#	Return converted format
			return (alst_vec, vlst_vec)
	end

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
	@doc raw"""
	**Description**  
	Simulates SIR/SIRS diffusion matching SAS exactly with preserved column order
	for RNG synchronization and cumulative peer effects.

	**Usage**  
	`sirdif(alst, vlst, infectedp, inf_r, rec_t, maxtime, p_symp, b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp; transmission_method=:weighted, immunity_duration=nothing)`

	**Arguments**
	- `alst`: Adjacency list (Matrix or Vector{Vector})
	- `vlst`: Edge weights aligned to `alst`
	- `infectedp`: Initial infected node IDs
	- `inf_r`: Base transmission probability
	- `rec_t`: Days to recovery (constant)
	- `maxtime`: Simulation horizon (days)
	- `p_symp`: Probability infected is symptomatic
	- `b_int, b_close, b_cxn_peers, b_cxn_total, b_cxn_symp, b_cls_x_smp`: Logit coefficients
	- `transmission_method`: `:weighted` or `:simple`
	- `immunity_duration`: Days of immunity after recovery (SIRS mode)

	**Details**
	- Preserves exact column order from adjacency matrix for RNG alignment
	- `nbrsinf` is cumulative and persistent across entire simulation
	- Infected processed in FIFO order by infection time
	- Immediate within-timestep updates when transmission occurs
	- **SAS-style termination**: when infections drop to zero, a single final row is written at `time=maxtime`.

	**Value**
	Dict with `infection_log` (time series), `total_time`, and `final_state`

	**See Also**
	`sim_prep`, `replicate_sas_simulation`
	""" sirdif

#   Exporting Objects
    export sim_prep,
		   sirdif

# 	Bring in the CLI submodule
	include("CLI.jl")   # defines module diffusion_sim.CLI

# 	Executable entrypoint for PackageCompiler:
"""
	Executable entrypoint. Returns 0 on success, non-zero on error.
"""
	function julia_main()::Cint
		try
			# Delegate to CLI with ARGS when the binary starts
			CLI.cli_main(ARGS)
			return 0
		catch err
			# Print a concise error + backtrace; non-zero exit for shell/R
			@error "Fatal error" exception=(err, catch_backtrace())
			return 1
		end
	end

end # module diffustion_sim

#Batch Simulation Tests: Seeds, UT-1.4, Replay, Seed-File Mode, UT-1.8, and SAS Network
#Jonathan H. Morgan, Ph.D.
#26 October 2026
#Revised 26 October 2026: outputs kept in sbi_tests/output; the SAS-transformed
#network option is tested. Seed components: weak components and eligible seeds on
#a hand-built network, and a batch on net_58 with seeds in the largest component.
#Run with several threads so that UT-1.8 compares serial and threaded batches:
#	julia -t 8 sbi_test.jl

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

################
#   PACKAGES   #
################

#	Loading Packages
	using JSON3
	using Random
	using Statistics
	using diffusion_sim
	using diffusion_sim.richards_model: TargetRegion

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Batch Tests: report a dictionary of checks
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
				println("  $(rpad(name, 46)) $(passed ? "PASS" : "FAIL")")
			end
			all_passed = all(values(checks))
			println("RESULT: $(all_passed ? "ALL CHECKS PASSED" : "SOME CHECKS FAILED")")

		#	Assembling Result
			return all_passed
	end

#	Helper Function for the Batch Tests: check that a call throws a given error type
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

#	Helper Function for the Batch Tests: read a JSONL file into parsed records
	function read_records(out_file::String)
		"""
		Args:
			out_file::String: JSONL output of run_sbi_batch
		Returns:
			NamedTuple: (lines::Vector{String}, records::Vector)
		"""

		#	Read and parse each line
			lines = readlines(out_file)
			records = [JSON3.read(line) for line in lines]

		#	Assembling Result
			return (lines = lines, records = records)
	end

#	Seeds: SplitMix64, Run Seeds, Eligible Nodes, and Seed Node Draws
	function test_seed_functions(network_data::NamedTuple)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
		Returns:
			Bool: true if every check passed
		"""

		#	SplitMix64 against its published first output for state 0
			checks = Dict{String, Bool}()
			checks["splitmix64(0) = 0xe220a8397b1dcdaf"] = diffusion_sim.splitmix64(UInt64(0)) == 0xe220a8397b1dcdaf

		#	Run seeds: deterministic, below 2^53, and distinct
			s1 = run_seeds(42, 7, 1)
			checks["run seeds are deterministic"] = s1 == run_seeds(42, 7, 1)
			seeds = [run_seeds(42, d, r) for d in 1:1000 for r in 1:3]
			run_vals = [s.run_seed for s in seeds]
			checks["run seeds lie in [0, 2^53)"] = all(0 .<= run_vals .< 2^53) && all(s -> 0 <= s.node_seed < 2^53, seeds)
			checks["3000 run seeds are distinct"] = length(unique(run_vals)) == 3000
			checks["run and node seeds differ"] = all(s -> s.run_seed != s.node_seed, seeds)
			checks["master seed changes the run seed"] = run_seeds(43, 7, 1).run_seed != s1.run_seed

		#	Eligible nodes have a tie to another node (any component)
			alst = network_data.alst
			eligible = eligible_seed_nodes(alst; component=:any)
			row_of = Dict(alst[i, 1] => i for i in 1:size(alst, 1))
			isolated = count(i -> !any(>(0), alst[i, 2:end]), 1:size(alst, 1))
			checks["every eligible node has a tie"] = all(id -> any(nbr -> nbr > 0 && nbr != id, alst[row_of[id], 2:end]), eligible)
			checks["eligible plus isolated is every node"] = length(eligible) + isolated == size(alst, 1)
			println("  ($(length(eligible)) eligible nodes, $(isolated) isolated)")

		#	Seed node draws: distinct, eligible, and validated
			rng = Xoshiro(1)
			draws = [draw_seed_nodes(rng, eligible, 3) for _ in 1:1000]
			eligible_set = Set(eligible)
			checks["seed nodes distinct within a draw"] = all(d -> length(unique(d)) == 3, draws)
			checks["seed nodes are eligible"] = all(d -> all(in(eligible_set), d), draws)
			checks["too many seed nodes rejected"] = throws_error(() -> draw_seed_nodes(rng, eligible, length(eligible) + 1), ArgumentError)

		#	Assembling Result
			return report_checks("Seeds: SplitMix64, Run Seeds, Eligible Nodes, Draws", checks)
	end

#	UT-1.4: sbi_batch End-to-End
	function test_ut14(graphml_path::String, params_file::String, work_dir::String)
		"""
		Args:
			graphml_path::String: network file
			params_file::String: the 20-draw test batch
			work_dir::String: directory for the outputs
		Returns:
			NamedTuple: (passed::Bool, out_file::String)
		Notes:
			The region is permissive (K_max = 1), so every established run is
			in the target. Criteria (roadmap UT-1.4): (a) 20 lines and a
			manifest; (b) every block and field present; (c) in_target null
			for non-starters, failed fits, and extinct runs under :exclude,
			and Boolean for established runs; (d) at least one run in the
			target; (e) draw_id values unique and in batch order; (f) no NaN
			or Inf in the file. Also checks that theta matches the input
			exactly and that the manifest records the settings.
		"""

		#	Run the batch
			out_file = joinpath(work_dir, "ut14.jsonl")
			res = run_sbi_batch(graphml_path, params_file, out_file; village="WeakCore1_3.2_2",
								region=TargetRegion(K_max=1.0), master_seed=42)
			out = read_records(out_file)
			input = JSON3.read(read(params_file, String))
			checks = Dict{String, Bool}()

		#	(a) Line count and manifest
			checks["a: 20 lines"] = length(out.lines) == 20 && res.n_runs == 20
			checks["a: manifest written"] = isfile(res.manifest_file)

		#	(b) Blocks and fields
			summary_fields = (:K, :r_growth, :t_i, :a, :frac_at_inflection, :fit_rmse, :peak_incidence,
							  :peak_prevalence, :peak_infected, :peak_day, :infected_at_surge, :max_slope,
							  :final_size, :duration, :extinguished, :t_10, :t_half, :t_90, :C0, :fit_ok,
							  :extinct, :viable, :status, :class)
			target_fields = (:in_target, :peak_ok, :surge_ok, :slope_ok, :prevalence_ok, :final_size_ok)
			meta_fields = (:village, :draw_id, :replicate, :run_seed, :seed_nodes, :n_nodes, :schema_version)
			fields_ok = all(out.records) do rec
				all(haskey(rec, b) for b in (:theta, :summary, :target, :meta)) &&
				all(haskey(rec[:summary], f) for f in summary_fields) &&
				all(haskey(rec[:target], f) for f in target_fields) &&
				all(haskey(rec[:meta], f) for f in meta_fields)
			end
			checks["b: every block and field present"] = fields_ok
			checks["b: schema version 3"] = all(rec -> rec[:meta][:schema_version] == 3, out.records)

		#	(c) in_target by class
			class_ok = all(out.records) do rec
				cls = rec[:summary][:class]
				it = rec[:target][:in_target]
				cls == "established" ? (it isa Bool) : (it === nothing)
			end
			checks["c: in_target null unless established"] = class_ok

		#	(d) At least one run in the target
			checks["d: at least one in_target = true"] = any(rec -> rec[:target][:in_target] === true, out.records)

		#	(e) draw_id unique and in batch order
			ids = [rec[:meta][:draw_id] for rec in out.records]
			checks["e: draw_id unique, batch order"] = ids == [d[:draw_id] for d in input]

		#	(f) No NaN or Inf in the file
			text = read(out_file, String)
			checks["f: no NaN or Inf"] = !occursin("NaN", text) && !occursin("Inf", text)

		#	theta equals the input exactly (imm null stays null)
			theta_fields = (:inf_r, :b_int, :b_close, :b_cxn_peers, :b_cxn_total, :b_cxn_symp,
							:b_cls_x_smp, :p_symp, :rec_scale, :rec_shape, :imm_scale, :imm_shape)
			theta_ok = all(zip(out.records, input)) do (rec, d)
				all(isequal(rec[:theta][f], d[f]) || rec[:theta][f] == d[f] for f in theta_fields)
			end
			checks["theta matches the input"] = theta_ok

		#	Manifest records the settings
			m = JSON3.read(read(res.manifest_file, String))
			checks["manifest records the settings"] = m[:master_seed] == 42 && m[:n_runs] == 20 &&
													  m[:maxtime] == 200 && m[:extinct_policy] == "exclude" &&
													  m[:target][:K_max] == 1.0 && m[:target][:I_peak_max] === nothing

		#	Report the classes
			println("  Classes: ", join(["$(k)=$(v)" for (k, v) in sort(collect(res.class_counts))], ", "))

		#	Assembling Result
			return (passed = report_checks("UT-1.4: sbi_batch End-to-End (20 draws)", checks), out_file = out_file)
	end

#	Replay: Stored Seeds Reproduce Records Exactly
	function test_replay(network_data::NamedTuple, params_file::String, out_file::String; n_replay::Int=5)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep (the same network as the batch)
			params_file::String: the batch's parameter file
			out_file::String: the batch output to replay from
			n_replay::Int: records to replay (default = 5)
		Returns:
			Bool: true if every replayed record is identical to the stored line
		Notes:
			Each replay calls simulate_and_summarize with the stored seed
			nodes and Xoshiro(run_seed), rebuilds the record, and compares
			the JSON text with the stored line (UT-1.8 criterion b).
		"""

		#	Inputs
			out = read_records(out_file)
			draws = diffusion_sim.read_draws(params_file)
			region = TargetRegion(K_max=1.0)
			checks = Dict{String, Bool}()

		#	Replay the first records
			for k in 1:min(n_replay, length(out.records))
				meta = out.records[k][:meta]
				seed_nodes = Int.(collect(meta[:seed_nodes]))
				summary = simulate_and_summarize(network_data.alst, network_data.vlst, seed_nodes, draws[k].theta,
												 200, region; rng=Xoshiro(meta[:run_seed]))
				rebuilt = (village = meta[:village], draw_id = meta[:draw_id], replicate = meta[:replicate],
						   run_seed = meta[:run_seed], seed_nodes = seed_nodes, n_nodes = meta[:n_nodes],
						   schema_version = meta[:schema_version])
				line = JSON3.write(diffusion_sim.build_record(draws[k].theta, summary, rebuilt))
				checks["draw $(meta[:draw_id]) replays exactly"] = line == out.lines[k]
			end

		#	Assembling Result
			return report_checks("Replay: Stored Seeds Reproduce Records", checks)
	end

#	Seed-File Mode Reproduces a Random-Mode Batch
	function test_seed_file_mode(graphml_path::String, params_file::String, out_file::String, work_dir::String)
		"""
		Args:
			graphml_path::String: network file
			params_file::String: the batch's parameter file
			out_file::String: a random-mode batch output
			work_dir::String: directory for the outputs
		Returns:
			Bool: true if every check passed
		Notes:
			write_seed_file extracts each run's seed and seed nodes; rerunning
			in testing mode must reproduce every line. Invalid seed files
			must be rejected.
		"""

		#	Write the seed file and rerun in file mode
			seed_file = joinpath(work_dir, "ut14_seeds.json")
			n_written = write_seed_file(out_file, seed_file)
			file_out = joinpath(work_dir, "ut14_file_mode.jsonl")
			run_sbi_batch(graphml_path, params_file, file_out; village="WeakCore1_3.2_2",
						  region=TargetRegion(K_max=1.0), seed_mode=:file, seed_file=seed_file)
			checks = Dict{String, Bool}()
			checks["seed file has 20 entries"] = n_written == 20
			checks["file mode reproduces every line"] = readlines(file_out) == readlines(out_file)

		#	A seed file missing a run is rejected
			short_file = joinpath(work_dir, "short_seeds.json")
			entries = JSON3.read(read(seed_file, String))
			open(short_file, "w") do io
				JSON3.write(io, collect(entries)[1:end-1])
			end
			checks["missing run rejected"] = throws_error(() -> run_sbi_batch(graphml_path, params_file,
														  joinpath(work_dir, "x.jsonl"); village="x",
														  seed_mode=:file, seed_file=short_file), ArgumentError)

		#	Settings errors are rejected before any run
			checks["random mode without seed rejected"] = throws_error(() -> run_sbi_batch(graphml_path, params_file,
																	   joinpath(work_dir, "x.jsonl"); village="x"), ArgumentError)
			checks["unknown policy rejected"] = throws_error(() -> run_sbi_batch(graphml_path, params_file,
															 joinpath(work_dir, "x.jsonl"); village="x", master_seed=1,
															 extinct_policy=:maybe), ArgumentError)

		#	Assembling Result
			return report_checks("Seed-File Mode and Settings Errors", checks)
	end

#	UT-1.8: Parallel Reproducibility
	function test_ut18(graphml_path::String, params_file::String, work_dir::String; replicates::Int=10)
		"""
		Args:
			graphml_path::String: network file
			params_file::String: the 20-draw test batch
			work_dir::String: directory for the outputs
			replicates::Int: replicates per draw (default = 10, so 200 runs)
		Returns:
			Bool: true if every check passed
		Notes:
			Runs the same batch serially and across all of Julia's threads.
			The two files must be identical byte for byte (UT-1.8 criterion
			a). With one thread the comparison is trivial, so the test warns
			and the check is reported as not exercised. Replicates must carry
			numbers 1 to R and distinct run seeds. Timing is informative.
		"""

		#	Serial and threaded batches
			serial_out = joinpath(work_dir, "ut18_serial.jsonl")
			thread_out = joinpath(work_dir, "ut18_threaded.jsonl")
			t_serial = @elapsed run_sbi_batch(graphml_path, params_file, serial_out; village="WeakCore1_3.2_2",
											  region=TargetRegion(K_max=1.0), master_seed=2026,
											  replicates=replicates, threaded=false)
			t_thread = @elapsed run_sbi_batch(graphml_path, params_file, thread_out; village="WeakCore1_3.2_2",
											  region=TargetRegion(K_max=1.0), master_seed=2026,
											  replicates=replicates, threaded=true)
			checks = Dict{String, Bool}()
			n_threads = Threads.nthreads()

		#	(a) Identical outputs
			identical = read(serial_out) == read(thread_out)
			if n_threads > 1
				checks["a: serial and $(n_threads)-thread files identical"] = identical
			else
				println("  WARNING: one thread available; start Julia with -t 8 to exercise UT-1.8")
				checks["a: identical (not exercised, 1 thread)"] = identical
			end

		#	Replicates: numbers 1 to R and distinct run seeds
			out = read_records(serial_out)
			reps = [rec[:meta][:replicate] for rec in out.records]
			seeds = [rec[:meta][:run_seed] for rec in out.records]
			checks["$(20 * replicates) runs written"] = length(out.records) == 20 * replicates
			checks["replicates numbered 1 to $(replicates)"] = reps == repeat(collect(1:replicates), 20)
			checks["run seeds distinct"] = length(unique(seeds)) == length(seeds)

		#	Report timing
			n_runs = 20 * replicates
			println("  Serial:   $(round(t_serial, digits=2)) s ($(round(n_runs / t_serial, digits=1)) runs/s)")
			println("  Threaded: $(round(t_thread, digits=2)) s ($(round(n_runs / t_thread, digits=1)) runs/s) " *
					"on $(n_threads) threads, speedup $(round(t_serial / t_thread, digits=2))")

		#	Assembling Result
			return report_checks("UT-1.8: Parallel Reproducibility ($(n_runs) runs)", checks)
	end

#	SAS-Transformed Network Option
	function test_sas_network(graphml_path::String, params_file::String, work_dir::String)
		"""
		Args:
			graphml_path::String: network file
			params_file::String: the 20-draw test batch
			work_dir::String: directory for the outputs
		Returns:
			Bool: true if every check passed
		Notes:
			Runs the batch with sas_transformation = true. The manifest must
			record the choice, n_nodes must match the transformed network,
			every seed node must be eligible in it, and a record must replay
			exactly when the network is loaded the same way. The run classes
			are reported for comparison with the raw network.
		"""

		#	Run the batch on the transformed network
			out_file = joinpath(work_dir, "ut14_sas.jsonl")
			res = run_sbi_batch(graphml_path, params_file, out_file; village="WeakCore1_3.2_2",
								sas_transformation=true, region=TargetRegion(K_max=1.0), master_seed=42)
			out = read_records(out_file)
			net_sas = sim_prep(graphml_path; sas_transformation=true)
			eligible = Set(eligible_seed_nodes(net_sas.alst))
			m = JSON3.read(read(res.manifest_file, String))
			checks = Dict{String, Bool}()

		#	Manifest, network size, and seeds
			checks["manifest records sas_transformation"] = m[:sas_transformation] === true
			checks["n_nodes matches the transformed network"] = all(rec -> rec[:meta][:n_nodes] == size(net_sas.alst, 1), out.records)
			checks["seed nodes eligible in the transformed network"] = all(rec -> all(in(eligible), rec[:meta][:seed_nodes]), out.records)

		#	Replay the first record on the transformed network
			draws = diffusion_sim.read_draws(params_file)
			meta = out.records[1][:meta]
			seed_nodes = Int.(collect(meta[:seed_nodes]))
			summary = simulate_and_summarize(net_sas.alst, net_sas.vlst, seed_nodes, draws[1].theta, 200,
											 TargetRegion(K_max=1.0); rng=Xoshiro(meta[:run_seed]))
			rebuilt = (village = meta[:village], draw_id = meta[:draw_id], replicate = meta[:replicate],
					   run_seed = meta[:run_seed], seed_nodes = seed_nodes, n_nodes = meta[:n_nodes],
					   schema_version = meta[:schema_version])
			line = JSON3.write(diffusion_sim.build_record(draws[1].theta, summary, rebuilt))
			checks["record replays on the transformed network"] = line == out.lines[1]

		#	Report the classes
			println("  Classes (SAS network): ", join(["$(k)=$(v)" for (k, v) in sort(collect(res.class_counts))], ", "))

		#	Assembling Result
			return report_checks("SAS-Transformed Network Option", checks)
	end

#	Seed Components: Weak Components, Eligible Seeds, and a Village Batch
	function test_seed_components(graphml_path::String, village_path::String, params_file::String, work_dir::String)
		"""
		Args:
			graphml_path::String: the connected test network (WeakCore1_3.2_2)
			village_path::String: a synthetic village with many components (net_58)
			params_file::String: the 20-draw test batch
			work_dir::String: directory for the outputs
		Returns:
			Bool: true if every check passed
		Notes:
			The hand-built network has IDs 11 to 18: {11, 12, 13, 14, 15} joined
			by ties in both directions and the one-directional tie 14 -> 11, with
			15 named by 11 but naming nobody; {16, 17}; and 18, whose only tie is
			to itself. On the connected test network the two rules must give the
			same eligible nodes in the same order, so earlier batches there are
			unchanged by the new default.
		"""

		#	Hand-built network: components and eligible nodes
			alst = [11 12 15; 12 13 0; 13 11 0; 14 11 0; 15 0 0; 16 17 0; 17 16 0; 18 18 0]
			comps = diffusion_sim.weak_components(alst)
			checks = Dict{String, Bool}()
			checks["toy: 3 weak components"] = length(comps.sizes) == 3
			checks["toy: largest has 5 nodes, row 1's set"] = comps.sizes[comps.largest] == 5 && comps.label[1] == comps.largest
			checks["toy: one-directional tie joins 14"] = comps.label[4] == comps.largest
			checks["toy: self-tie leaves 18 alone"] = comps.sizes[comps.label[8]] == 1
			checks["toy: largest eligible = [11, 12, 13, 14]"] = eligible_seed_nodes(alst) == [11, 12, 13, 14]
			checks["toy: any eligible = [11, 12, 13, 14, 16, 17]"] = eligible_seed_nodes(alst; component=:any) == [11, 12, 13, 14, 16, 17]
			checks["toy: unknown component rejected"] = throws_error(() -> eligible_seed_nodes(alst; component=:all), ArgumentError)

		#	Equal-sized components: the one with the earliest row is the largest
			tied = [21 22 0; 22 21 0; 23 24 0; 24 23 0]
			checks["toy: size tie goes to the earliest row"] = eligible_seed_nodes(tied) == [21, 22]

		#	Connected test network: both rules agree, raw and SAS-transformed
			for sas in (false, true)
				net = sim_prep(graphml_path; sas_transformation=sas)
				checks["connected network (sas=$(sas)): largest == any"] =
					eligible_seed_nodes(net.alst) == eligible_seed_nodes(net.alst; component=:any)
			end

		#	Village batch with seeds in the largest component
			out_file = joinpath(work_dir, "net_58_largest.jsonl")
			res = run_sbi_batch(village_path, params_file, out_file; village="net_58", sas_transformation=true,
								region=TargetRegion(K_max=1.0), master_seed=42)
			out = read_records(out_file)
			m = JSON3.read(read(res.manifest_file, String))
			net = sim_prep(village_path; sas_transformation=true)
			vc = diffusion_sim.weak_components(net.alst)
			row_of = Dict(net.alst[i, 1] => i for i in 1:size(net.alst, 1))
			in_largest = id -> vc.label[row_of[id]] == vc.largest
			checks["net_58: manifest seed_component = largest"] = m[:seed_component] == "largest"
			checks["net_58: manifest component counts"] = m[:n_components] == length(vc.sizes) &&
														  m[:largest_component] == vc.sizes[vc.largest]
			checks["net_58: manifest n_eligible_seeds"] = m[:n_eligible_seeds] == length(eligible_seed_nodes(net.alst))
			checks["net_58: every seed in the largest component"] = all(rec -> all(in_largest, rec[:meta][:seed_nodes]), out.records)

		#	The same batch with seeds from any component, for comparison
			res_any = run_sbi_batch(village_path, params_file, joinpath(work_dir, "net_58_any.jsonl"); village="net_58",
									sas_transformation=true, seed_component=:any, region=TargetRegion(K_max=1.0),
									master_seed=42)
			println("  net_58: $(m[:n_components]) components, largest $(m[:largest_component]) of $(m[:n_nodes]) nodes, ",
					"$(m[:n_eligible_seeds]) eligible seeds")
			println("  Classes, largest: ", join(["$(k)=$(v)" for (k, v) in sort(collect(res.class_counts))], ", "))
			println("  Classes, any:     ", join(["$(k)=$(v)" for (k, v) in sort(collect(res_any.class_counts))], ", "))

		#	Assembling Result
			return report_checks("Seed Components: Weak Components and Eligible Seeds", checks)
	end

#############
#   TESTS   #
#############

#	Paths (EC2 test data) and a scratch directory for outputs
	graphml_path = "/workspace/ergm_data/sim_test_data/WeakCore1_3.2_2.graphml"
	params_file = "/workspace/ergm-based-diffusion-simulations/sbi_tests/params_batch_ut14.json"
	village_path = "/workspace/ergm_data/sim_test_data/synthetic_village_m_graphml/net_58.graphml"
	work_dir = mkpath("/workspace/ergm-based-diffusion-simulations/sbi_tests/output")
	network_data = sim_prep(graphml_path)

#	Warm-up (compilation) so that UT-1.8 timing excludes it
	run_sbi_batch(graphml_path, params_file, joinpath(work_dir, "warmup.jsonl"); village="warmup", master_seed=1)

#	Run Each Test
	seeds_passed = test_seed_functions(network_data)
	ut14 = test_ut14(graphml_path, params_file, work_dir)
	replay_passed = test_replay(network_data, params_file, ut14.out_file)
	file_mode_passed = test_seed_file_mode(graphml_path, params_file, ut14.out_file, work_dir)
	ut18_passed = test_ut18(graphml_path, params_file, work_dir)
	sas_passed = test_sas_network(graphml_path, params_file, work_dir)
	components_passed = test_seed_components(graphml_path, village_path, params_file, work_dir)

#	Summary
	println("=" ^ 60)
	println("Batch Simulation Test Summary")
	println("=" ^ 60)
	println("  Seed functions:            $(seeds_passed ? "PASS" : "FAIL")")
	println("  UT-1.4 sbi_batch:          $(ut14.passed ? "PASS" : "FAIL")")
	println("  Replay from stored seeds:  $(replay_passed ? "PASS" : "FAIL")")
	println("  Seed-file mode:            $(file_mode_passed ? "PASS" : "FAIL")")
	println("  UT-1.8 parallel:           $(ut18_passed ? "PASS" : "FAIL")")
	println("  SAS network option:        $(sas_passed ? "PASS" : "FAIL")")
	println("  Seed components:           $(components_passed ? "PASS" : "FAIL")")
	println("  Outputs in: $(work_dir)")

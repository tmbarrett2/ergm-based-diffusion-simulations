#sirdif Equivalence Test: Changes to sirdif's Internals Against the Frozen Reference
#Jonathan H. Morgan, Ph.D.
#26 October 2026
#Protocol: run this script once with the current diffusion_sim.jl, which confirms
#that reference_sirdif_v2.jl matches it, then replace diffusion_sim.jl with the
#changed version and run it again. Both runs must pass. Re-run simulator_tests.jl
#(UT-1.7) after the replacement.

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

################
#   PACKAGES   #
################

#	Loading Packages
	using DataFrames
	using Random
	using Statistics
	using diffusion_sim

#   Frozen Copy of sirdif with Weibull Durations and Waning Immunity
	include("/workspace/ergm-based-diffusion-simulations/sbi_tests/reference_sirdif_v2.jl")

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Equivalence Tests: report a dictionary of checks
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

#	Helper Function for the Equivalence Tests: build a small network in sim_prep format
	function build_test_network(neighbors::Vector{Vector{Int}})
		"""
		Args:
			neighbors::Vector{Vector{Int}}: neighbors[i] lists the IDs that node i names (IDs are 1:n)
		Returns:
			NamedTuple: (alst::Matrix{Int}, vlst::Matrix{Float64}) in sim_prep format
		Notes:
			Unlike the builder in simulator_tests.jl, ties need not be listed
			in both directions, and a node may name itself, so the edge
			cases of one-directional ties and self-ties can be tested.
			All edge weights are 1.0.
		"""

		#	Size the matrices
			n = length(neighbors)
			max_degree = max(1, maximum(length.(neighbors)))
			alst = zeros(Int, n, max_degree + 1)
			vlst = zeros(Float64, n, max_degree + 1)

		#	Fill IDs, neighbors, and weights
			for i in 1:n
				alst[i, 1] = i
				vlst[i, 1] = Float64(i)
				for (col, j) in enumerate(neighbors[i])
					alst[i, col + 1] = j
					vlst[i, col + 1] = 1.0
				end
			end

		#	Assembling Result
			return (alst = alst, vlst = vlst)
	end

#	Helper Function for the Equivalence Tests: compare two sirdif results exactly
	function identical_runs(new_res::Dict, ref_res::Dict)
		"""
		Args:
			new_res::Dict: result from diffusion_sim.sirdif
			ref_res::Dict: result from ReferenceSirdifV2.sirdif
		Returns:
			Bool: true if the logs and final states are identical
		Notes:
			Every column of infection_log, including prop_cum_infected, and
			every column of final_state must match element by element.
			isequal treats two NaN values as equal, which == does not.
		"""

		#	Logs and final states
			log_match = isequal(new_res["infection_log"], ref_res["infection_log"])
			state_match = new_res["final_state"] == ref_res["final_state"]

		#	Assembling Result
			return log_match && state_match
	end

#	Equivalence on the Test Network
	function test_equivalence_network(network_data::NamedTuple; n_seeds::Int=100, maxtime::Int=200)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
			n_seeds::Int: runs per setting (default = 100)
			maxtime::Int: simulation horizon (default = 200)
		Returns:
			Bool: true if every run in every setting is identical
		Notes:
			Four settings with the R-setting coefficients: Weibull recovery
			with Weibull immunity; short, variable recovery with short
			immunity, which produces the most losses of immunity; fixed
			recovery with fixed immunity; and Weibull recovery with permanent
			immunity. Run s of each setting gives both simulators their own
			Xoshiro(s), and draws its seed node from a separate generator.
		"""

		#	Network and eligible seed nodes
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)
			ids = alst[:, 1]
			eligible = findall(i -> alst[i, 2] > 0, 1:n)

		#	Settings: (name, recovery, immunity)
			settings = [("Weibull recovery, Weibull immunity", WeibullDuration(10.0, 2.0), WeibullDuration(30.0, 2.0)),
						("short recovery, short immunity", WeibullDuration(3.0, 0.8), WeibullDuration(8.0, 1.5)),
						("fixed recovery, fixed immunity", FixedDuration(14), FixedDuration(20)),
						("Weibull recovery, permanent immunity", WeibullDuration(10.0, 2.0), nothing)]
			checks = Dict{String, Bool}()

		#	Compare the simulators run by run
			for (name, recovery, immunity) in settings
				n_match = 0
				for s in 1:n_seeds
					seed_id = ids[eligible[rand(Xoshiro(1_000_000 + s), 1:length(eligible))]]
					ref_res = ReferenceSirdifV2.sirdif(alst, vlst, [seed_id], 0.33, recovery, maxtime, 0.5,
													   -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
													   immunity=immunity, rng=Xoshiro(s))
					new_res = sirdif(alst, vlst, [seed_id], 0.33, recovery, maxtime, 0.5,
									 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1;
									 immunity=immunity, rng=Xoshiro(s))
					n_match += identical_runs(new_res, ref_res)
				end
				checks["$(name) ($(n_match)/$(n_seeds))"] = n_match == n_seeds
			end

		#	Assembling Result
			return report_checks("Equivalence on the Test Network ($(n) nodes)", checks)
	end

#	Equivalence on Edge Cases
	function test_equivalence_edge_cases(; n_seeds::Int=200, maxtime::Int=60)
		"""
		Args:
			n_seeds::Int: runs per seed node (default = 200)
			maxtime::Int: simulation horizon (default = 60)
		Returns:
			Bool: true if every run is identical
		Notes:
			An 8-node network with a self-tie (node 3 names itself), a
			one-directional chain (4 names 5, 5 names 6, but not the reverse),
			and an isolated node (8). High activation, high transmission,
			short recovery, and two-day immunity make infections, recoveries,
			and losses of immunity frequent. Runs are seeded at nodes 1, 3,
			and 4 in turn.
		"""

		#	Network with the edge cases
			neighbors = [[2, 3], [1, 3, 4], [1, 2, 3], [5], [6], [5, 7], [6], Int[]]
			net = build_test_network(neighbors)
			checks = Dict{String, Bool}()

		#	Compare the simulators run by run
			for seed_node in [1, 3, 4]
				n_match = 0
				for s in 1:n_seeds
					ref_res = ReferenceSirdifV2.sirdif(net.alst, net.vlst, [seed_node], 0.5, WeibullDuration(2.0, 1.5),
													   maxtime, 0.5, 2.0, 0.5, 0.0, 0.0, 0.0, 0.0;
													   immunity=FixedDuration(2), rng=Xoshiro(s))
					new_res = sirdif(net.alst, net.vlst, [seed_node], 0.5, WeibullDuration(2.0, 1.5),
									 maxtime, 0.5, 2.0, 0.5, 0.0, 0.0, 0.0, 0.0;
									 immunity=FixedDuration(2), rng=Xoshiro(s))
					n_match += identical_runs(new_res, ref_res)
				end
				checks["seeded at node $(seed_node) ($(n_match)/$(n_seeds))"] = n_match == n_seeds
			end

		#	Assembling Result
			return report_checks("Equivalence on Edge Cases (self-tie, one-directional ties, isolate)", checks)
	end

#	Timing and Allocation: Reference Against Current sirdif
	function test_timing(network_data::NamedTuple; n_runs::Int=100, maxtime::Int=200)
		"""
		Args:
			network_data::NamedTuple: output from sim_prep
			n_runs::Int: runs per simulator (default = 100)
			maxtime::Int: simulation horizon (default = 200)
		Returns:
			NamedTuple: (ref_ms, new_ms, ref_mib, new_mib) per run
		Notes:
			Informative, not a pass criterion. Uses the setting of the
			26 October timing check (Weibull recovery and immunity, one seed
			node, R-setting coefficients). Both simulators are compiled by a
			warm-up run before timing.
		"""

		#	Network and eligible seed nodes
			alst = network_data.alst
			vlst = network_data.vlst
			n = size(alst, 1)
			ids = alst[:, 1]
			eligible = findall(i -> alst[i, 2] > 0, 1:n)
			seed_ids = [ids[eligible[rand(Xoshiro(2026 + r), 1:length(eligible))]] for r in 1:n_runs]
			rec = WeibullDuration(10.0, 2.0)
			imm = WeibullDuration(30.0, 2.0)

		#	Warm-up
			ReferenceSirdifV2.sirdif(alst, vlst, [seed_ids[1]], 0.33, rec, maxtime, 0.5,
									 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1; immunity=imm, rng=Xoshiro(1))
			sirdif(alst, vlst, [seed_ids[1]], 0.33, rec, maxtime, 0.5,
				   -0.1, 1.0, -0.5, -3.5, -1.5, -0.1; immunity=imm, rng=Xoshiro(1))

		#	Reference simulator
			ref = @timed for r in 1:n_runs
				ReferenceSirdifV2.sirdif(alst, vlst, [seed_ids[r]], 0.33, rec, maxtime, 0.5,
										 -0.1, 1.0, -0.5, -3.5, -1.5, -0.1; immunity=imm, rng=Xoshiro(2026 + r))
			end

		#	Current simulator
			cur = @timed for r in 1:n_runs
				sirdif(alst, vlst, [seed_ids[r]], 0.33, rec, maxtime, 0.5,
					   -0.1, 1.0, -0.5, -3.5, -1.5, -0.1; immunity=imm, rng=Xoshiro(2026 + r))
			end

		#	Report per run
			result = (ref_ms = 1000 * ref.time / n_runs, new_ms = 1000 * cur.time / n_runs,
					  ref_mib = ref.bytes / n_runs / 2^20, new_mib = cur.bytes / n_runs / 2^20)
			println("=" ^ 60)
			println("Timing and Allocation per Run ($(n_runs) runs, informative)")
			println("=" ^ 60)
			println("  reference sirdif: $(round(result.ref_ms, digits=2)) ms, $(round(result.ref_mib, digits=2)) MiB")
			println("  current sirdif:   $(round(result.new_ms, digits=2)) ms, $(round(result.new_mib, digits=2)) MiB")

		#	Assembling Result
			return result
	end

#############
#   TESTS   #
#############

#	Network Path (EC2 test data)
	graphml_path = "/workspace/ergm_data/sim_test_data/WeakCore1_3.2_2.graphml"
	network_data = sim_prep(graphml_path)

#	Run Each Test
	network_passed = test_equivalence_network(network_data)
	edge_passed = test_equivalence_edge_cases()
	timing = test_timing(network_data)

#	Summary
	println("=" ^ 60)
	println("sirdif Equivalence Test Summary")
	println("=" ^ 60)
	println("  Test network (four settings):  $(network_passed ? "PASS" : "FAIL")")
	println("  Edge cases:                    $(edge_passed ? "PASS" : "FAIL")")
	println("  Time per run, reference:       $(round(timing.ref_ms, digits=2)) ms")
	println("  Time per run, current:         $(round(timing.new_ms, digits=2)) ms")

#   SocFeedback Simulation Precompilation Script
#   Jonathan H. Morgan, Ph.D.
#   Revised 26 October 2026: runs a small sbi_batch through the CLI so that the
#   compiled binary includes the batch, Richards fit, and JSON code paths.

#   Load Diffusion Sim
    using diffusion_sim

#   Load Common Dependencies
    using CSV, DataFrames, EzXML, StatsBase, JSON3, ArgParse, ProgressMeter, Statistics, Random

#   Verify that the precompile script has loaded
    @info "precompile_exec.jl loaded"

#	Small Network and Parameter Draws for the Precompile Batch
#	A 30-node ring with ties to the second neighbor, written as undirected GraphML,
#	and four draws chosen so that runs spread and the Richards fit is exercised
	work_dir = mktempdir()
	graphml_file = joinpath(work_dir, "precompile_ring.graphml")
	params_file = joinpath(work_dir, "precompile_params.json")
	out_file = joinpath(work_dir, "precompile_results.jsonl")

	io = open(graphml_file, "w")
	try
		println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
		println(io, "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\">")
		println(io, "  <graph id=\"G\" edgedefault=\"undirected\">")
		for i in 1:30
			println(io, "    <node id=\"n$(i)\"/>")
		end
		for i in 1:30
			println(io, "    <edge source=\"n$(i)\" target=\"n$(mod1(i + 1, 30))\"/>")
			println(io, "    <edge source=\"n$(i)\" target=\"n$(mod1(i + 2, 30))\"/>")
		end
		println(io, "  </graph>")
		println(io, "</graphml>")
	finally
		close(io)
	end

	draws = [(draw_id = k, inf_r = 0.5, b_int = 2.0, b_close = 1.0, b_cxn_peers = -0.5, b_cxn_total = -1.0,
			  b_cxn_symp = -1.0, b_cls_x_smp = -0.1, p_symp = 0.5, rec_scale = 5.0, rec_shape = 2.0,
			  imm_scale = k == 4 ? nothing : 30.0, imm_shape = k == 4 ? nothing : 2.0) for k in 1:4]
	io = open(params_file, "w")
	try
		JSON3.write(io, draws)
	finally
		close(io)
	end

#	Run the Batch Through the CLI Entry Point
	diffusion_sim.CLI.cli_main(["sbi_batch", "--graphml", graphml_file, "--village", "precompile",
								"--params", params_file, "--out", out_file, "--max-time", "60",
								"--seed", "1", "--target-K-max", "1.0"])
	@info "precompile sbi_batch finished: $(countlines(out_file)) records"

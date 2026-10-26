#Batch Simulation Tests: Seeds, UT-1.4, Replay, Seed-File Mode, and UT-1.8
#Jonathan H. Morgan, Ph.D.
#26 October 2026
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

#################
#   FUNCTIONS   #
#################

#	Helper Function for the Comparison: flatten a record into dotted keys
	function flatten_record!(flat::Dict{String, Any}, prefix::String, x)
		if x isa JSON3.Object
			for (k, v) in pairs(x)
				flatten_record!(flat, prefix == "" ? String(k) : "$(prefix).$(k)", v)
			end
		else
			flat[prefix] = x isa JSON3.Array ? collect(x) : x
		end
		return flat
	end

###########################
#   COMPARE JSONL FILES   #
###########################

#	Compare every line and field
	a = readlines(ARGS[1]); b = readlines(ARGS[2])
	println("Lines: $(length(a)) vs $(length(b))")
	worst = Dict{String, Float64}(); exact_diffs = Dict{String, Int}()
	for i in 1:min(length(a), length(b))
		fa = flatten_record!(Dict{String, Any}(), "", JSON3.read(a[i]))
		fb = flatten_record!(Dict{String, Any}(), "", JSON3.read(b[i]))
		for k in union(keys(fa), keys(fb))
			va = get(fa, k, missing); vb = get(fb, k, missing)
			isequal(va, vb) && continue
			if va isa AbstractFloat && vb isa AbstractFloat
				worst[k] = max(get(worst, k, 0.0), abs(va - vb) / max(abs(va), abs(vb)))
			else
				exact_diffs[k] = get(exact_diffs, k, 0) + 1
			end
		end
	end

#	Report
	println("Numeric fields that differ (largest relative difference):")
	for (k, v) in sort(collect(worst), by = last, rev = true)
		println("  $(rpad(k, 32)) $(v)")
	end
	println("Non-numeric or exact fields that differ (number of runs):")
	isempty(exact_diffs) ? println("  none") : for (k, v) in sort(collect(exact_diffs)); println("  $(rpad(k, 32)) $(v)"); end
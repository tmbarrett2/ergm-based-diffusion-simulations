#SIR Social Feedback Simulation
#Jonathan H. Morgan, Ph.D., James Moody Ph.D., & Tyler Barret
#27 August 2025

#   Pulling-In diffustion_sim & Activating Local Environment
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate("/workspace/ergm-based-diffusion-simulations/julia_env")
    Pkg.status()

#   Setting the DISPLAY
    ENV["DISPLAY"] = ":100"

################
#   PACKAGES   #
################
using CSV
using DataFrames
using Statistics
using diffustion_sim

######################
#   TEST FUNCTIONS   #
######################

#	Quick test function for specific parameter combinations
	function quick_param_test(b_int::Float64, b_cxn_symp::Float64)
		"""
		Args:
			b_int::Float64: baseline activity parameter to test
			b_cxn_symp::Float64: symptomatic reduction parameter to test
		Returns:
			DataFrame: subset of results for quick inspection
		Notes:
			Runs sim_parm_space with reduced parameter space for quick testing.
		"""
		
		#	Run with reduced parameter space
			results = sim_parm_space(
				b_int,
				b_cxn_symp,
				-2.0,   # moderate global effect
				-0.3,   # moderate peer effect
				0.5,    # moderate edge boost
				-0.05,  # small interaction term
				2,      # just 2 edge weights
				3,      # just 3 max peers
				500
			)
			
		#	Return summary
			println("Quick test with b_int=", b_int, ", b_cxn_symp=", b_cxn_symp)
			println("Generated ", nrow(results), " rows")
			println("Probability range: [", 
				round(minimum(results.prob_act), digits=4), ", ",
				round(maximum(results.prob_act), digits=4), "]")
			
			return results
	end

######################
#   FUNCTION TESTS   #
######################

#   arithmetic_mode 

#   Test 1: Simple Vector (PASSED)
    test1 = [1, 2, 3, 2, 2, 4, 5]
    result1 = arithmetic_mode(test1)
    println("Test 1 - Input: $test1")
    println("Mode: $result1")
    println("Expected: 2")
    println()

#   Test 2: Array with missing values (keeping missing)
    test2 = [1, 2, missing, 2, 3, missing, missing]
    result2 = arithmetic_mode(test2, na_rm=false)
    println("Test 2 - Input with missing (na_rm=false): $test2")
    println("Mode: $result2")
    println("Expected: missing (most frequent)")
    println()

#   Test 3: Array with missing values (removing missing)
    test3 = [1, 2, missing, 2, 3, missing, 2]
    result3 = arithmetic_mode(test3, na_rm=true)
    println("Test 3 - Input with missing (na_rm=true): $test3")
    println("Mode: $result3")
    println("Expected: 2")
    println()

#   Test 4: Float array
    test4 = [1.5, 2.3, 1.5, 3.7, 1.5, 2.3]
    result4 = arithmetic_mode(test4)
    println("Test 4 - Float input: $test4")
    println("Mode: $result4")
    println("Expected: 1.5")
    println()

#   Test 5: String array
    test5 = ["apple", "banana", "apple", "cherry", "apple", "banana"]
    result5 = arithmetic_mode(test5)
    println("Test 5 - String input: $test5")
    println("Mode: $result5")
    println("Expected: apple")

#   color_assignment

#   Test 1: Simple test with values across the range
    test1 = [0.02, 0.07, 0.12, 0.18, 0.23, 0.28, 0.33, 0.38, 0.43, 0.48, 0.52, 0.58]
    result1 = color_assignment(test1)
 

    println("Test 1 - Input values: $test1")
    println("Colors (optimized): $result1")
    println()

#   Test 2: Edge cases
    test2 = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55]
    result2 = color_assignment(test2)
    println("Test 2 - Edge values: $test2")
    println("Colors: $result2")
    println()

#   Test 3: Random values
    test3 = rand(10) * 0.6  # Random values between 0 and 0.6
    result3 = color_assignment(test3)
    println("Test 3 - Random values: $(round.(test3, digits=3))")
    println("Colors: $result3")
    println()

#   Test 4: Single value
    test4 = [0.27]
    result4 = color_assignment(test4)
    println("Test 4 - Single value: $test4")
    println("Color (optimized): $result4")
    println("Expected: navajowhite (for range [0.25, 0.30))")
    println()

#   Parameter Space Exploration

#   Full test with output
    test_results = test_sim_parm_space(verbose=true, export_csv=false)

#   Plot results
    results = sim_parm_space(-0.1, -1.5, -3.5, -0.5, 1.0, -0.1, 3, 5, 500)
    plot_parameter_space(results, edgemax=3)

#   Load network from GraphML file
	network_data = sim_prep("/workspace/data/sim_test_data/WeakCore1_3.2_2.graphml")
	alst = network_data.alst
	vlst = network_data.vlst
	
#   Check network structure
	n_nodes = size(alst, 1)
	max_degree = size(alst, 2) - 1
	println("Network has $n_nodes nodes with max degree $max_degree")

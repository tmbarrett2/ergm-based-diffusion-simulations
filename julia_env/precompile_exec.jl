#   SocFeedback Simulation Precompilation Script
#   Jonathan H. Morgan, Ph.D.

#   Load Diffusion Sim
    using diffusion_sim

#   Load Common Dependencies
    using CSV, DataFrames, EzXML, StatsBase, JSON3, ArgParse, ProgressMeter, Statistics, Random

#   Verify that the precompile script has loaded
    @info "precompile_exec.jl loaded"

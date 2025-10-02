#   Building the SocFeedback_Simulator's Executable
#   Jonathan H. Morgan, Ph.D.
#   2 October 2025

#   Activate Package
    cd("/workspace/ergm-based-diffusion-simulations")
    using Pkg
    Pkg.activate(@__DIR__)
    Pkg.status()
    Pkg.resolve()
    Pkg.instantiate()

    ENV["JULIA_DEBUG"] = "PackageCompiler"
    ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0"

    using PackageCompiler

#   Check the Entry Point
    @info "Verifying module & entrypoint…"
    import diffusion_sim
    @assert isdefined(diffusion_sim, :julia_main) "diffusion_sim.julia_main is not defined"

#   Specify Parameters/Constants
    const PKGDIR = normpath(@__DIR__)
    const APPDIR = joinpath(@__DIR__, "build", "diffusion_sim_app")
    const PRECOMPILE_EXEC = joinpath(@__DIR__, "precompile_exec.jl")
    @info "Checking precompile script at $PRECOMPILE_EXEC"
    @assert isfile(PRECOMPILE_EXEC) "Missing precompile_exec.jl at $PRECOMPILE_EXEC"

#   Build the App
    create_app(
        PKGDIR,
        APPDIR;
        executables = ["diffusion_sim" => "julia_main"],  # <-- point to the function
        precompile_execution_file = PRECOMPILE_EXEC,
        incremental = false,
        force = true,
        filter_stdlibs = false,
        include_transitive_dependencies = true,
        include_preferences = true,
        cpu_target = "generic",
        sysimage_build_args = `--compile=all`
    )

#   Notices
    @info "Done. Run it like:"
    @info joinpath(APPDIR, "bin", "diffusion_sim") * " --help"

    println("\n✅ Built app at: ", APPDIR)
    println("Binary:\n  ", joinpath(APPDIR, "bin", "diffusion_sim"))
    println("\nTry:")
    println("  ", joinpath(APPDIR, "bin", "diffusion_sim"), " --help")

# Julia Compilation Command # julia --project -e 'include("build_app.jl")' 
# ./julia_env/build/diffusion_sim_app/bin/diffusion_sim --help

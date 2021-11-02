import TurbulenceConvection
using Profile
using Test

tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "Cases.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "parameter_set.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "main.jl"))
import .NameList

TurbulenceConvection.initialize_io(sim::Simulation1d) = nothing
TurbulenceConvection.io(sim::Simulation1d) = nothing

function run_main(; time_run = false)
    case_name = "Bomex"
    println("Running $case_name...")
    namelist = NameList.default_namelist(case_name)
    namelist["time_stepping"]["t_max"] = namelist["time_stepping"]["t_max"] / 10
    namelist["stats_io"]["frequency"] = 10.0e10
    namelist["meta"]["uuid"] = "01"
    main(namelist; time_run = time_run)
end

# run_main(; time_run = true)
run_main() # run first to compile
# @profile run_main()

Profile.print()

import Plots
# import FlameGraphs
import ProfileView

ProfileView.@profview run_main()
# g = FlameGraphs.flamegraph(C=true)
# img = FlameGraphs.flamepixels(g);
# Plots.plot(img)

Profile.clear_malloc_data()
nothing


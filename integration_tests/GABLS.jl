if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test
using Random

# Make deterministic:
Random.seed!(1234)

include(joinpath("utils", "Cases.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .Cases
using .NameList
using .ParamList

include(joinpath("utils", "main.jl"))

best_mse = OrderedDict()

best_mse["updraft_area"] = 1.8227678974387011e+01
best_mse["updraft_w"] = 1.2435239287489708e+01

best_mse["updraft_thetal"] = 4.6453137166954967e+01
best_mse["v_mean"] = 1.0939547963673784e+01
best_mse["u_mean"] = 5.7078250617016124e+02
best_mse["tke_mean"] = 7.7344960671102756e+00

ds_pycles = Dataset(joinpath(PyCLES_output_dataset_path, "Gabls.nc"), "r")

@testset "GABLS" begin
    println("Running GABLS...")
    namelist = NameList.GABLS(default_namelist("GABLS"))
    paramlist = ParamList.GABLS(default_paramlist("GABLS"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        compute_mse(
            ds,
            ds_pycles,
            "GABLS",
            best_mse,
            dirname(ds_filename);
            plot_comparison=true
        )
    end

    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end


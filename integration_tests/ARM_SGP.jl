if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_paramlist.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList
using .ParamList

best_mse = OrderedDict()
best_mse["qt_mean"] = 6.7939652328911526e+00
best_mse["updraft_area"] = 2.6728907094911216e+02
best_mse["updraft_w"] = 3.3246977024115265e-01
best_mse["updraft_qt"] = 6.0684798609715308e+01
best_mse["updraft_thetal"] = 6.9056681132480222e+01
best_mse["u_mean"] = 8.7994360629093748e+01
best_mse["tke_mean"] = 4.1453532207769657e+00

@testset "ARM_SGP" begin
    println("Running ARM_SGP...")
    namelist = NameList.ARM_SGP(default_namelist("ARM_SGP"))
    paramlist = ParamList.ARM_SGP(default_paramlist("ARM_SGP"))
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist, paramlist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "ARM_SGP.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "ARM_SGP.nc"), "r") do ds_scampy
                compute_mse(
                    "ARM_SGP",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    nothing
end

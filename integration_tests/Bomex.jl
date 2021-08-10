if !("." in LOAD_PATH) # for easier local testing
    push!(LOAD_PATH, ".")
end
import TurbulenceConvection
using TurbulenceConvection
using Test

include(joinpath("utils", "main.jl"))
include(joinpath("utils", "generate_namelist.jl"))
include(joinpath("utils", "compute_mse.jl"))
using .NameList

best_mse = OrderedDict()
best_mse["qt_mean"] = 9.5377282321874049e-02
best_mse["updraft_area"] = 7.2522420468798543e+02
best_mse["updraft_w"] = 2.6282677670021386e+01
best_mse["updraft_qt"] = 4.0341139540222635e+00
best_mse["updraft_thetal"] = 2.1548164483602367e+01
best_mse["v_mean"] = 6.5247595142753596e+01
best_mse["u_mean"] = 5.3292989796462209e+01
best_mse["tke_mean"] = 3.8891129116688681e+01
best_mse["temperature_mean"] = 3.8547176978568602e-05
best_mse["ql_mean"] = 3.8134281154981489e+00
best_mse["thetal_mean"] = 3.9162071364386570e-05
best_mse["Hvar_mean"] = 7.5510623226760771e+01
best_mse["QTvar_mean"] = 2.2962565404185188e+01

@testset "Bomex" begin
    println("Running Bomex...")
    namelist = default_namelist("Bomex")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Bomex.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Bomex.nc"), "r") do ds_scampy
                compute_mse(
                    "Bomex",
                    best_mse,
                    joinpath(dirname(ds_filename), "comparison");
                    ds_turb_conv = ds,
                    ds_scampy = ds_scampy,
                    ds_pycles = ds_pycles,
                    plot_comparison = true,
                    t_start = 4 * 3600,
                    t_stop = 6 * 3600,
                )
            end
        end
    end

    for k in keys(best_mse)
        test_mse(computed_mse, best_mse, k)
    end
    nothing
end

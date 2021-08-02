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
best_mse["qt_mean"] = 8.2239735839064965e-02
best_mse["updraft_area"] = 7.0282567418077099e+02
best_mse["updraft_w"] = 8.7997106745562263e+01
best_mse["updraft_qt"] = 5.9534668549222083e+00
best_mse["updraft_thetal"] = 2.3059746455183749e+01
best_mse["v_mean"] = 1.2344255867930939e+02
best_mse["u_mean"] = 5.3486981021128493e+01
best_mse["tke_mean"] = 3.2038984505851317e+01
best_mse["temperature_mean"] = 3.1019627325766795e-05
best_mse["ql_mean"] = 1.5630034030829634e+01
best_mse["thetal_mean"] = 3.1601763199735002e-05
best_mse["Hvar_mean"] = 3.4265143495072330e+01
best_mse["QTvar_mean"] = 1.4543323453279942e+01

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
                )
            end
        end
    end

    test_mse(computed_mse, best_mse, "qt_mean")
    test_mse(computed_mse, best_mse, "updraft_area")
    test_mse(computed_mse, best_mse, "updraft_w")
    test_mse(computed_mse, best_mse, "updraft_qt")
    test_mse(computed_mse, best_mse, "updraft_thetal")
    test_mse(computed_mse, best_mse, "v_mean")
    test_mse(computed_mse, best_mse, "u_mean")
    test_mse(computed_mse, best_mse, "tke_mean")
    test_mse(computed_mse, best_mse, "temperature_mean")
    test_mse(computed_mse, best_mse, "ql_mean")
    test_mse(computed_mse, best_mse, "thetal_mean")
    nothing
end

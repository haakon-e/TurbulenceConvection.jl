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
best_mse["qt_mean"] = 5.4596292300803628e-01
best_mse["updraft_area"] = 1.6644272370900098e+03
best_mse["updraft_w"] = 5.5689731617925463e+02
best_mse["updraft_qt"] = 2.4446559332424545e+01
best_mse["updraft_thetal"] = 6.5601697606639760e+01
best_mse["v_mean"] = 1.0610811546766800e+02
best_mse["u_mean"] = 1.1387870745418634e+02
best_mse["tke_mean"] = 9.1707254861508807e+02
best_mse["temperature_mean"] = 2.4294131149581795e-04
best_mse["ql_mean"] = 1.3122452320829298e+04
best_mse["thetal_mean"] = 1.8634680078510443e-04
best_mse["Hvar_mean"] = 1.9466083203062441e+04
best_mse["QTvar_mean"] = 1.0274739422376710e+05

@testset "Rico" begin
    println("Running Rico...")
    namelist = default_namelist("Rico")
    namelist["meta"]["uuid"] = "01"
    ds_filename = @time main(namelist)

    computed_mse = Dataset(ds_filename, "r") do ds
        Dataset(joinpath(PyCLES_output_dataset_path, "Rico.nc"), "r") do ds_pycles
            Dataset(joinpath(SCAMPy_output_dataset_path, "Rico.nc"), "r") do ds_scampy
                compute_mse(
                    "Rico",
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
    test_mse(computed_mse, best_mse, "Hvar_mean")
    test_mse(computed_mse, best_mse, "QTvar_mean")
    nothing
end

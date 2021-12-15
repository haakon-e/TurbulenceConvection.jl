import TurbulenceConvection
import Plots
import NCDatasets
import CLIMAParameters
import ClimaCore
const CC = ClimaCore
tc_dir = dirname(dirname(pathof(TurbulenceConvection)))
include(joinpath(tc_dir, "integration_tests", "utils", "generate_namelist.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "Cases.jl"))
include(joinpath(tc_dir, "integration_tests", "utils", "parameter_set.jl"))
using .NameList
import .Cases
function export_ref_profile(case_name::String)
    TC = TurbulenceConvection
    namelist = default_namelist(case_name)
    param_set = create_parameter_set(namelist)
    namelist["meta"]["uuid"] = "01"

    FT = Float64
    grid = TC.Grid(FT(namelist["grid"]["dz"]), namelist["grid"]["nz"])
    Stats = TC.NetCDFIO_Stats(namelist, grid)
    case = Cases.get_case(namelist)
    ref_params = Cases.reference_params(case, grid, param_set, namelist)
    ref_state = TC.ReferenceState(grid, param_set, Stats; ref_params...)

    aux_vars(FT) = (; ref_state = (ρ0 = FT(0), α0 = FT(0), p0 = FT(0)))
    aux_cent_fields = TC.FieldFromNamedTuple(TC.center_space(grid), aux_vars(FT))
    aux_face_fields = TC.FieldFromNamedTuple(TC.face_space(grid), aux_vars(FT))

    parent(aux_face_fields.ref_state.p0) .= ref_state.p0
    parent(aux_face_fields.ref_state.ρ0) .= ref_state.rho0
    parent(aux_face_fields.ref_state.α0) .= ref_state.alpha0
    parent(aux_cent_fields.ref_state.p0) .= ref_state.p0_half
    parent(aux_cent_fields.ref_state.ρ0) .= ref_state.rho0_half
    parent(aux_cent_fields.ref_state.α0) .= ref_state.alpha0_half

    aux = CC.Fields.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    io_nt = TC.io_dictionary_ref_state((; aux))
    TC.initialize_io(io_nt, Stats)
    TC.io(io_nt, Stats)

    NCDatasets.Dataset(joinpath(Stats.path_plus_file), "r") do ds
        zc = ds.group["profiles"]["zc"][:]
        zf = ds.group["profiles"]["zf"][:]
        ρ0_c = ds.group["reference"]["ρ0_c"][:]
        p0_c = ds.group["reference"]["p0_c"][:]
        α0_c = ds.group["reference"]["α0_c"][:]
        ρ0_f = ds.group["reference"]["ρ0_f"][:]
        p0_f = ds.group["reference"]["p0_f"][:]
        α0_f = ds.group["reference"]["α0_f"][:]

        p1 = Plots.plot(ρ0_c, zc ./ 1000; label = "centers")
        Plots.plot!(ρ0_f, zf ./ 1000; label = "faces")
        Plots.plot!(size = (1000, 400))
        Plots.plot!(margin = 5 * Plots.mm)
        Plots.xlabel!("ρ_0")
        Plots.ylabel!("z (km)")
        Plots.title!("ρ_0")

        p2 = Plots.plot(p0_c ./ 1000, zc ./ 1000; label = "centers")
        Plots.plot!(p0_f ./ 1000, zf ./ 1000; label = "faces")
        Plots.plot!(size = (1000, 400))
        Plots.plot!(margin = 5 * Plots.mm)
        Plots.xlabel!("p_0 (kPa)")
        Plots.ylabel!("z (km)")
        Plots.title!("p_0 (kPa)")

        p3 = Plots.plot(α0_c, zc ./ 1000; label = "centers")
        Plots.plot!(α0_f, zf ./ 1000; label = "faces")
        Plots.plot!(size = (1000, 400))
        Plots.plot!(margin = 5 * Plots.mm)
        Plots.xlabel!("α_0")
        Plots.ylabel!("z (km)")
        Plots.title!("α_0")
        Plots.plot(p1, p2, p3; layout = (1, 3))
        Plots.savefig("$case_name.svg")
    end
end
for case_name in ("Bomex", "life_cycle_Tan2018", "Soares", "Rico", "ARM_SGP", "DYCOMS_RF01", "GABLS", "SP", "DryBubble")
    export_ref_profile(case_name)
end;

# Note: temperatures in this case become extremely low.
CLIMAParameters.Planet.T_freeze(::EarthParameterSet) = 100.0
export_ref_profile("TRMM_LBA")
export_ref_profile("GATE_III")

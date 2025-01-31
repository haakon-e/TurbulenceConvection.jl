import TurbulenceConvection
const TC = TurbulenceConvection
import Thermodynamics
const TD = Thermodynamics

function initialize_edmf(edmf::TC.EDMF_PrognosticTKE, grid, state, Case, gm::TC.GridMeanVariables, TS::TC.TimeStepping)
    initialize_covariance(edmf, grid, state, gm, Case)
    up = edmf.UpdVar
    param_set = TC.parameter_set(gm)
    aux_tc = TC.center_aux_turbconv(state)
    prog_gm = TC.center_prog_grid_mean(state)
    p0_c = TC.center_ref_state(state).p0
    parent(aux_tc.prandtl_nvec) .= edmf.prandtl_number
    @inbounds for k in TC.real_center_indices(grid)
        ts = TC.thermo_state_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        aux_tc.θ_virt[k] = TD.virtual_pottemp(ts)
    end
    TC.update_surface(Case, grid, state, gm, TS, param_set)
    TC.compute_updraft_surface_bc(edmf, grid, state, Case)
    if Case.casename == "DryBubble"
        initialize_updrafts_DryBubble(edmf, grid, state, up, gm)
    else
        initialize_updrafts(edmf, grid, state, up, gm)
    end
    TC.set_edmf_surface_bc(edmf, grid, state, up, Case.Sur)
    return
end

function initialize_covariance(edmf::TC.EDMF_PrognosticTKE, grid, state, gm, Case)

    kc_surf = TC.kc_surface(grid)
    en = edmf.EnvVar
    aux_gm = TC.center_aux_grid_mean(state)
    prog_en = TC.center_prog_environment(state)
    aux_en = TC.center_aux_environment(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    aux_bulk = TC.center_aux_bulk(state)
    ae = 1 .- aux_bulk.area # area of environment

    aux_en.tke .= aux_gm.tke
    prog_en.ρatke .= aux_en.tke .* ρ0_c .* ae

    TC.get_GMV_CoVar(edmf, grid, state, :tke, :w)
    aux_gm.Hvar .= aux_gm.Hvar[kc_surf] .* aux_gm.tke
    aux_gm.QTvar .= aux_gm.QTvar[kc_surf] .* aux_gm.tke
    aux_gm.HQTcov .= aux_gm.HQTcov[kc_surf] .* aux_gm.tke

    prog_en.ρaHvar .= aux_gm.Hvar .* ρ0_c .* ae
    prog_en.ρaQTvar .= aux_gm.QTvar .* ρ0_c .* ae
    prog_en.ρaHQTcov .= aux_gm.HQTcov .* ρ0_c .* ae
    return
end

function initialize_updrafts(edmf, grid, state, up::TC.UpdraftVariables, gm::TC.GridMeanVariables)
    kc_surf = TC.kc_surface(grid)
    aux_up = TC.center_aux_updrafts(state)
    prog_gm = TC.center_prog_grid_mean(state)
    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_up = TC.center_prog_updrafts(state)
    prog_up_f = TC.face_prog_updrafts(state)
    ρ0_c = TC.center_ref_state(state).ρ0
    @inbounds for i in 1:(up.n_updrafts)
        @inbounds for k in TC.real_face_indices(grid)
            aux_up_f[i].w[k] = 0
            prog_up_f[i].ρaw[k] = 0
        end

        @inbounds for k in TC.real_center_indices(grid)
            aux_up[i].buoy[k] = 0
            # Simple treatment for now, revise when multiple updraft closures
            # become more well defined
            aux_up[i].area[k] = 0
            aux_up[i].q_tot[k] = prog_gm.q_tot[k]
            aux_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
            aux_up[i].q_liq[k] = aux_gm.q_liq[k]
            aux_up[i].q_ice[k] = aux_gm.q_ice[k]
            aux_up[i].T[k] = aux_gm.T[k]
            prog_up[i].ρarea[k] = 0
            prog_up[i].ρaq_tot[k] = 0
            prog_up[i].ρaθ_liq_ice[k] = 0
        end

        aux_up[i].area[kc_surf] = edmf.area_surface_bc[i]
        prog_up[i].ρarea[kc_surf] = ρ0_c[kc_surf] * aux_up[i].area[kc_surf]
    end
    return
end

function initialize_updrafts_DryBubble(edmf, grid, state, up::TC.UpdraftVariables, gm::TC.GridMeanVariables)
    dz = grid.Δz

    # criterion 2: b>1e-4
    #! format: off
    z_in = [
          75.,  125.,  175.,  225.,  275.,  325.,  375.,  425.,  475.,
         525.,  575.,  625.,  675.,  725.,  775.,  825.,  875.,  925.,
         975., 1025., 1075., 1125., 1175., 1225., 1275., 1325., 1375.,
        1425., 1475., 1525., 1575., 1625., 1675., 1725., 1775., 1825.,
        1875., 1925., 1975., 2025., 2075., 2125., 2175., 2225., 2275.,
        2325., 2375., 2425., 2475., 2525., 2575., 2625., 2675., 2725.,
        2775., 2825., 2875., 2925., 2975., 3025., 3075., 3125., 3175.,
        3225., 3275., 3325., 3375., 3425., 3475., 3525., 3575., 3625.,
        3675., 3725., 3775., 3825., 3875., 3925.]

    θ_liq_in = [
        299.9882, 299.996 , 300.0063, 300.0205, 300.04  , 300.0594,
        300.0848, 300.1131, 300.1438, 300.1766, 300.2198, 300.2567,
        300.2946, 300.3452, 300.3849, 300.4245, 300.4791, 300.5182,
        300.574 , 300.6305, 300.6668, 300.7222, 300.7771, 300.8074,
        300.8591, 300.9092, 300.9574, 300.9758, 301.0182, 301.0579,
        301.0944, 301.1276, 301.1572, 301.1515, 301.1729, 301.1902,
        301.2033, 301.2122, 301.2167, 301.2169, 301.2127, 301.2041,
        301.1913, 301.1743, 301.1533, 301.1593, 301.1299, 301.097 ,
        301.0606, 301.0212, 300.9788, 300.9607, 300.9125, 300.8625,
        300.8108, 300.7806, 300.7256, 300.6701, 300.6338, 300.5772,
        300.5212, 300.482 , 300.4272, 300.3875, 300.3354, 300.2968,
        300.2587, 300.2216, 300.1782, 300.1452, 300.1143, 300.0859,
        300.0603, 300.0408, 300.0211, 300.0067, 299.9963, 299.9884    ]

    Area_in = [
        0.04 , 0.055, 0.07 , 0.08 , 0.085, 0.095, 0.1  , 0.105, 0.11 ,
        0.115, 0.115, 0.12 , 0.125, 0.125, 0.13 , 0.135, 0.135, 0.14 ,
        0.14 , 0.14 , 0.145, 0.145, 0.145, 0.15 , 0.15 , 0.15 , 0.15 ,
        0.155, 0.155, 0.155, 0.155, 0.155, 0.155, 0.16 , 0.16 , 0.16 ,
        0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 , 0.16 ,
        0.155, 0.155, 0.155, 0.155, 0.155, 0.155, 0.15 , 0.15 , 0.15 ,
        0.15 , 0.145, 0.145, 0.145, 0.14 , 0.14 , 0.14 , 0.135, 0.135,
        0.13 , 0.13 , 0.125, 0.12 , 0.115, 0.115, 0.11 , 0.105, 0.1  ,
        0.095, 0.085, 0.08 , 0.07 , 0.055, 0.04    ]

    W_in = [
        0.017 , 0.0266, 0.0344, 0.0417, 0.0495, 0.0546, 0.061 , 0.0668,
        0.0721, 0.0768, 0.0849, 0.0887, 0.092 , 0.0996, 0.1019, 0.1037,
        0.1106, 0.1114, 0.1179, 0.1243, 0.1238, 0.1297, 0.1355, 0.1335,
        0.1387, 0.1437, 0.1485, 0.1448, 0.1489, 0.1527, 0.1564, 0.1597,
        0.1628, 0.1565, 0.1588, 0.1609, 0.1626, 0.1641, 0.1652, 0.166 ,
        0.1665, 0.1667, 0.1666, 0.1662, 0.1655, 0.1736, 0.1722, 0.1706,
        0.1686, 0.1664, 0.1639, 0.1698, 0.1667, 0.1634, 0.1599, 0.1641,
        0.1601, 0.1559, 0.1589, 0.1543, 0.1496, 0.1514, 0.1464, 0.1475,
        0.1422, 0.1425, 0.1424, 0.1419, 0.1361, 0.135 , 0.1335, 0.1316,
        0.1294, 0.1302, 0.1271, 0.1264, 0.1269, 0.1256    ]

    T_in = [
        299.2557, 298.775 , 298.2969, 297.8227, 297.3536, 296.8843,
        296.421 , 295.9603, 295.502 , 295.0456, 294.5994, 294.1468,
        293.6951, 293.2556, 292.8054, 292.3549, 291.9188, 291.4677,
        291.0325, 290.5978, 290.1434, 289.7073, 289.2706, 288.81  ,
        288.3698, 287.928 , 287.4842, 287.0118, 286.5622, 286.1099,
        285.6544, 285.1957, 284.7335, 284.2379, 283.7677, 283.2937,
        282.8157, 282.3337, 281.8476, 281.3574, 280.8631, 280.3649,
        279.8626, 279.3565, 278.8467, 278.362 , 277.8447, 277.3241,
        276.8006, 276.2742, 275.7454, 275.2388, 274.705 , 274.1694,
        273.6327, 273.1155, 272.576 , 272.0363, 271.514 , 270.9736,
        270.4339, 269.9094, 269.3711, 268.8465, 268.311 , 267.7877,
        267.2649, 266.7432, 266.2159, 265.698 , 265.1821, 264.6685,
        264.1574, 263.6518, 263.1461, 262.6451, 262.1476, 261.6524]
    #! format: on

    aux_up = TC.center_aux_updrafts(state)
    aux_up_f = TC.face_aux_updrafts(state)
    aux_gm = TC.center_aux_grid_mean(state)
    prog_gm = TC.center_prog_grid_mean(state)
    prog_up = TC.center_prog_updrafts(state)
    ρ_0_c = TC.center_ref_state(state).ρ0
    ρ_0_f = TC.face_ref_state(state).ρ0
    Area_in = TC.pyinterp(grid.zc, z_in, Area_in)
    θ_liq_in = TC.pyinterp(grid.zc, z_in, θ_liq_in)
    T_in = TC.pyinterp(grid.zc, z_in, T_in)
    @inbounds for i in 1:(up.n_updrafts)
        @inbounds for k in TC.real_face_indices(grid)
            if minimum(z_in) <= grid.zf[k] <= maximum(z_in)
                aux_up_f[i].w[k] = 0.0
            end
        end

        @inbounds for k in TC.real_center_indices(grid)
            if minimum(z_in) <= grid.zc[k] <= maximum(z_in)
                aux_up[i].area[k] = Area_in[k]
                aux_up[i].θ_liq_ice[k] = θ_liq_in[k]
                aux_up[i].q_tot[k] = 0.0
                aux_up[i].q_liq[k] = 0.0
                aux_up[i].q_ice[k] = 0.0

                # for now temperature is provided as diagnostics from LES
                aux_up[i].T[k] = T_in[k]
                prog_up[i].ρarea[k] = ρ_0_c[k] * aux_up[i].area[k]
                prog_up[i].ρaθ_liq_ice[k] = prog_up[i].ρarea[k] * aux_up[i].θ_liq_ice[k]
                prog_up[i].ρaq_tot[k] = prog_up[i].ρarea[k] * aux_up[i].q_tot[k]
            else
                aux_up[i].area[k] = 0.0
                aux_up[i].θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
                aux_up[i].T[k] = aux_gm.T[k]
                prog_up[i].ρarea[k] = 0.0
                prog_up[i].ρaθ_liq_ice[k] = 0.0
                prog_up[i].ρaq_tot[k] = 0.0
            end
        end
    end
    return
end

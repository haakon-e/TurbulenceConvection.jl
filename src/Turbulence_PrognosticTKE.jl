
function initialize(edmf::EDMF_PrognosticTKE, grid, state, Case::CasesBase, gm::GridMeanVariables, TS::TimeStepping)
    initialize_covariance(edmf, grid, state, gm, Case)
    if Case.casename == "DryBubble"
        initialize_DryBubble(edmf, grid, state, edmf.UpdVar, gm)
    else
        initialize(edmf, grid, state, edmf.UpdVar, gm)
    end
    return
end

# Initialize the IO pertaining to this class
function initialize_io(edmf::EDMF_PrognosticTKE, Stats::NetCDFIO_Stats)

    initialize_io(edmf.UpdVar, Stats)
    initialize_io(edmf.EnvVar, Stats)
    initialize_io(edmf.Rain, Stats)

    add_profile(Stats, "eddy_viscosity")
    add_profile(Stats, "eddy_diffusivity")
    add_profile(Stats, "entrainment_sc")
    add_profile(Stats, "detrainment_sc")
    add_profile(Stats, "nh_pressure")
    add_profile(Stats, "nh_pressure_adv")
    add_profile(Stats, "nh_pressure_drag")
    add_profile(Stats, "nh_pressure_b")
    add_profile(Stats, "asp_ratio")
    add_profile(Stats, "horiz_K_eddy")
    add_profile(Stats, "sorting_function")
    add_profile(Stats, "b_mix")
    add_ts(Stats, "rd")
    add_profile(Stats, "turbulent_entrainment")
    add_profile(Stats, "turbulent_entrainment_full")
    add_profile(Stats, "turbulent_entrainment_W")
    add_profile(Stats, "turbulent_entrainment_H")
    add_profile(Stats, "turbulent_entrainment_QT")
    add_profile(Stats, "massflux")
    add_profile(Stats, "massflux_h")
    add_profile(Stats, "massflux_qt")
    add_profile(Stats, "massflux_tendency_h")
    add_profile(Stats, "massflux_tendency_qt")
    add_profile(Stats, "diffusive_flux_h")
    add_profile(Stats, "diffusive_flux_u")
    add_profile(Stats, "diffusive_flux_v")
    add_profile(Stats, "diffusive_flux_qt")
    add_profile(Stats, "diffusive_tendency_h")
    add_profile(Stats, "diffusive_tendency_qt")
    add_profile(Stats, "total_flux_h")
    add_profile(Stats, "total_flux_qt")
    add_profile(Stats, "mixing_length")
    add_profile(Stats, "updraft_qt_precip")
    add_profile(Stats, "updraft_thetal_precip")
    # Diff mixing lengths: Ignacio
    add_profile(Stats, "ed_length_scheme")
    add_profile(Stats, "mixing_length_ratio")
    add_profile(Stats, "entdet_balance_length")
    add_profile(Stats, "interdomain_tke_t")
    add_profile(Stats, "tke_buoy")
    add_profile(Stats, "tke_dissipation")
    add_profile(Stats, "tke_entr_gain")
    add_profile(Stats, "tke_detr_loss")
    add_profile(Stats, "tke_shear")
    add_profile(Stats, "tke_pressure")
    add_profile(Stats, "tke_interdomain")
    add_profile(Stats, "Hvar_dissipation")
    add_profile(Stats, "QTvar_dissipation")
    add_profile(Stats, "HQTcov_dissipation")
    add_profile(Stats, "Hvar_entr_gain")
    add_profile(Stats, "QTvar_entr_gain")
    add_profile(Stats, "Hvar_detr_loss")
    add_profile(Stats, "QTvar_detr_loss")
    add_profile(Stats, "HQTcov_detr_loss")
    add_profile(Stats, "HQTcov_entr_gain")
    add_profile(Stats, "Hvar_shear")
    add_profile(Stats, "QTvar_shear")
    add_profile(Stats, "HQTcov_shear")
    add_profile(Stats, "Hvar_rain")
    add_profile(Stats, "QTvar_rain")
    add_profile(Stats, "HQTcov_rain")
    add_profile(Stats, "Hvar_interdomain")
    add_profile(Stats, "QTvar_interdomain")
    add_profile(Stats, "HQTcov_interdomain")
    return
end

function io(edmf::EDMF_PrognosticTKE, grid, state, Stats::NetCDFIO_Stats, TS::TimeStepping, param_set)

    mean_nh_pressure = face_field(grid)
    mean_nh_pressure_adv = face_field(grid)
    mean_nh_pressure_drag = face_field(grid)
    mean_nh_pressure_b = face_field(grid)

    mean_asp_ratio = center_field(grid)
    mean_entr_sc = center_field(grid)
    mean_detr_sc = center_field(grid)
    massflux = center_field(grid)
    mean_frac_turb_entr = center_field(grid)
    mean_horiz_K_eddy = center_field(grid)
    mean_sorting_function = center_field(grid)
    mean_b_mix = center_field(grid)

    io(edmf.UpdVar, grid, state, Stats)
    io(edmf.EnvVar, grid, state, Stats)
    io(edmf.Rain, grid, state, Stats, edmf.UpdThermo, edmf.EnvThermo, TS)

    prog_up = center_prog_updrafts(state)
    aux_tc = center_aux_tc(state)
    a_up_bulk = aux_tc.bulk.area

    write_profile(Stats, "eddy_viscosity", diffusivity_m(edmf).values)
    write_profile(Stats, "eddy_diffusivity", diffusivity_h(edmf).values)
    write_ts(Stats, "rd", StatsBase.mean(edmf.pressure_plume_spacing))

    @inbounds for k in real_center_indices(grid)
        if a_up_bulk[k] > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                massflux[k] += interpf2c(edmf.m, grid, k, i)
                mean_entr_sc[k] += prog_up[i].area[k] * edmf.entr_sc[i, k] / a_up_bulk[k]
                mean_detr_sc[k] += prog_up[i].area[k] * edmf.detr_sc[i, k] / a_up_bulk[k]
                mean_asp_ratio[k] += prog_up[i].area[k] * edmf.asp_ratio[i, k] / a_up_bulk[k]
                mean_frac_turb_entr[k] += prog_up[i].area[k] * edmf.frac_turb_entr[i, k] / a_up_bulk[k]
                mean_horiz_K_eddy[k] += prog_up[i].area[k] * edmf.horiz_K_eddy[i, k] / a_up_bulk[k]
                mean_sorting_function[k] += prog_up[i].area[k] * edmf.sorting_function[i, k] / a_up_bulk[k]
                mean_b_mix[k] += prog_up[i].area[k] * edmf.b_mix[i, k] / a_up_bulk[k]
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        a_up_bulk_f =
            interpc2f(a_up_bulk, grid, k; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        if a_up_bulk_f > 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                a_up_f = interpc2f(
                    prog_up[i].area,
                    grid,
                    k;
                    bottom = SetValue(edmf.area_surface_bc[i]),
                    top = SetZeroGradient(),
                )
                mean_nh_pressure[k] += a_up_f * edmf.nh_pressure[i, k] / a_up_bulk_f
                mean_nh_pressure_b[k] += a_up_f * edmf.nh_pressure_b[i, k] / a_up_bulk_f
                mean_nh_pressure_adv[k] += a_up_f * edmf.nh_pressure_adv[i, k] / a_up_bulk_f
                mean_nh_pressure_drag[k] += a_up_f * edmf.nh_pressure_drag[i, k] / a_up_bulk_f
            end
        end
    end

    write_profile(Stats, "turbulent_entrainment", mean_frac_turb_entr)
    write_profile(Stats, "horiz_K_eddy", mean_horiz_K_eddy)
    write_profile(Stats, "entrainment_sc", mean_entr_sc)
    write_profile(Stats, "detrainment_sc", mean_detr_sc)
    write_profile(Stats, "nh_pressure", mean_nh_pressure)
    write_profile(Stats, "nh_pressure_adv", mean_nh_pressure_adv)
    write_profile(Stats, "nh_pressure_drag", mean_nh_pressure_drag)
    write_profile(Stats, "nh_pressure_b", mean_nh_pressure_b)
    write_profile(Stats, "asp_ratio", mean_asp_ratio)
    write_profile(Stats, "massflux", massflux)
    write_profile(Stats, "massflux_h", edmf.massflux_h)
    write_profile(Stats, "massflux_qt", edmf.massflux_qt)
    write_profile(Stats, "massflux_tendency_h", edmf.massflux_tendency_h)
    write_profile(Stats, "massflux_tendency_qt", edmf.massflux_tendency_qt)
    write_profile(Stats, "diffusive_flux_h", edmf.diffusive_flux_h)
    write_profile(Stats, "diffusive_flux_qt", edmf.diffusive_flux_qt)
    write_profile(Stats, "diffusive_flux_u", edmf.diffusive_flux_u)
    write_profile(Stats, "diffusive_flux_v", edmf.diffusive_flux_v)
    write_profile(Stats, "diffusive_tendency_h", edmf.diffusive_tendency_h)
    write_profile(Stats, "diffusive_tendency_qt", edmf.diffusive_tendency_qt)
    write_profile(Stats, "total_flux_h", edmf.massflux_h .+ edmf.diffusive_flux_h)
    write_profile(Stats, "total_flux_qt", edmf.massflux_h .+ edmf.diffusive_flux_qt)
    write_profile(Stats, "mixing_length", edmf.mixing_length)
    write_profile(Stats, "updraft_qt_precip", edmf.UpdThermo.qt_tendency_rain_formation_tot)
    write_profile(Stats, "updraft_thetal_precip", edmf.UpdThermo.θ_liq_ice_tendency_rain_formation_tot)

    #Different mixing lengths : Ignacio
    write_profile(Stats, "ed_length_scheme", edmf.mls)
    write_profile(Stats, "mixing_length_ratio", edmf.ml_ratio)
    write_profile(Stats, "entdet_balance_length", edmf.l_entdet)
    write_profile(Stats, "interdomain_tke_t", edmf.b)
    compute_covariance_dissipation(edmf, grid, state, edmf.EnvVar.TKE, param_set)
    write_profile(Stats, "tke_dissipation", edmf.EnvVar.TKE.dissipation)
    write_profile(Stats, "tke_entr_gain", edmf.EnvVar.TKE.entr_gain)
    compute_covariance_detr(edmf, grid, state, edmf.EnvVar.TKE)
    write_profile(Stats, "tke_detr_loss", edmf.EnvVar.TKE.detr_loss)
    write_profile(Stats, "tke_shear", edmf.EnvVar.TKE.shear)
    write_profile(Stats, "tke_buoy", edmf.EnvVar.TKE.buoy)
    write_profile(Stats, "tke_pressure", edmf.EnvVar.TKE.press)
    write_profile(Stats, "tke_interdomain", edmf.EnvVar.TKE.interdomain)

    compute_covariance_dissipation(edmf, grid, state, edmf.EnvVar.Hvar, param_set)
    write_profile(Stats, "Hvar_dissipation", edmf.EnvVar.Hvar.dissipation)
    compute_covariance_dissipation(edmf, grid, state, edmf.EnvVar.QTvar, param_set)
    write_profile(Stats, "QTvar_dissipation", edmf.EnvVar.QTvar.dissipation)
    compute_covariance_dissipation(edmf, grid, state, edmf.EnvVar.HQTcov, param_set)
    write_profile(Stats, "HQTcov_dissipation", edmf.EnvVar.HQTcov.dissipation)
    write_profile(Stats, "Hvar_entr_gain", edmf.EnvVar.Hvar.entr_gain)
    write_profile(Stats, "QTvar_entr_gain", edmf.EnvVar.QTvar.entr_gain)
    write_profile(Stats, "HQTcov_entr_gain", edmf.EnvVar.HQTcov.entr_gain)
    compute_covariance_detr(edmf, grid, state, edmf.EnvVar.Hvar)
    compute_covariance_detr(edmf, grid, state, edmf.EnvVar.QTvar)
    compute_covariance_detr(edmf, grid, state, edmf.EnvVar.HQTcov)
    write_profile(Stats, "Hvar_detr_loss", edmf.EnvVar.Hvar.detr_loss)
    write_profile(Stats, "QTvar_detr_loss", edmf.EnvVar.QTvar.detr_loss)
    write_profile(Stats, "HQTcov_detr_loss", edmf.EnvVar.HQTcov.detr_loss)
    write_profile(Stats, "Hvar_shear", edmf.EnvVar.Hvar.shear)
    write_profile(Stats, "QTvar_shear", edmf.EnvVar.QTvar.shear)
    write_profile(Stats, "HQTcov_shear", edmf.EnvVar.HQTcov.shear)
    write_profile(Stats, "Hvar_rain", edmf.EnvVar.Hvar.rain_src)
    write_profile(Stats, "QTvar_rain", edmf.EnvVar.QTvar.rain_src)
    write_profile(Stats, "HQTcov_rain", edmf.EnvVar.HQTcov.rain_src)
    write_profile(Stats, "Hvar_interdomain", edmf.EnvVar.Hvar.interdomain)
    write_profile(Stats, "QTvar_interdomain", edmf.EnvVar.QTvar.interdomain)
    write_profile(Stats, "HQTcov_interdomain", edmf.EnvVar.HQTcov.interdomain)
    return
end

#= These methods are to be overloaded by Cases.jl =#
function update_surface end
function update_forcing end
function update_radiation end

function update_cloud_frac(edmf::EDMF_PrognosticTKE, grid, state, GMV::GridMeanVariables)
    # update grid-mean cloud fraction and cloud cover
    prog_up = center_prog_updrafts(state)
    aux_tc = center_aux_tc(state)
    a_up_bulk = aux_tc.bulk.area
    @inbounds for k in real_center_indices(grid) # update grid-mean cloud fraction and cloud cover
        edmf.EnvVar.Area.values[k] = 1 - a_up_bulk[k]
        GMV.cloud_fraction.values[k] =
            edmf.EnvVar.Area.values[k] * edmf.EnvVar.cloud_fraction.values[k] +
            a_up_bulk[k] * edmf.UpdVar.cloud_fraction[k]
    end
    GMV.cloud_cover = min(edmf.EnvVar.cloud_cover + sum(edmf.UpdVar.cloud_cover), 1)
end

function compute_gm_tendencies!(edmf::EDMF_PrognosticTKE, grid, state, Case, gm, TS)
    gm.U.tendencies .= 0
    gm.V.tendencies .= 0
    gm.QT.tendencies .= 0
    gm.H.tendencies .= 0
    param_set = parameter_set(gm)
    ρ0_f = face_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    α0_c = center_ref_state(state).α0
    kf_surf = kf_surface(grid)
    kc_surf = kc_surface(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area # area of environment

    @inbounds for k in real_center_indices(grid)
        # Apply large-scale horizontal advection tendencies
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], gm.H.values[k], gm.QT.values[k])
        Π = TD.exner(ts)

        if Case.Fo.apply_coriolis
            gm.U.tendencies[k] -= Case.Fo.coriolis_param * (Case.Fo.vg[k] - gm.V.values[k])
            gm.V.tendencies[k] += Case.Fo.coriolis_param * (Case.Fo.ug[k] - gm.U.values[k])
        end
        if rad_type(Case.Rad) <: Union{RadiationDYCOMS_RF01, RadiationLES}
            gm.H.tendencies[k] += Case.Rad.dTdt[k] / Π
        end
        H_cut = ccut_downwind(gm.H.values, grid, k)
        q_tot_cut = ccut_downwind(gm.QT.values, grid, k)
        ∇H = c∇_downwind(H_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))
        ∇q_tot = c∇_downwind(q_tot_cut, grid, k; bottom = FreeBoundary(), top = SetGradient(0))

        if force_type(Case.Fo) <: ForcingDYCOMS_RF01
            gm.QT.tendencies[k] += Case.Fo.dqtdt[k]
            # Apply large-scale subsidence tendencies
            gm.H.tendencies[k] -= ∇H * Case.Fo.subsidence[k]
            gm.QT.tendencies[k] -= ∇q_tot * Case.Fo.subsidence[k]
        end

        if force_type(Case.Fo) <: ForcingStandard
            if Case.Fo.apply_subsidence
                gm.H.tendencies[k] -= ∇H * Case.Fo.subsidence[k]
                gm.QT.tendencies[k] -= ∇q_tot * Case.Fo.subsidence[k]
            end
            gm.H.tendencies[k] += Case.Fo.dTdt[k] / Π
            gm.QT.tendencies[k] += Case.Fo.dqtdt[k]
        end

        if force_type(Case.Fo) <: ForcingLES
            H_horz_adv = Case.Fo.dtdt_hadv[k] / Π
            H_nudge = Case.Fo.dtdt_nudge[k] / Π
            H_fluc = Case.Fo.dtdt_fluc[k] / Π

            gm_U_nudge_k = (Case.Fo.u_nudge[k] - gm.U.values[k]) / Case.Fo.nudge_tau
            gm_V_nudge_k = (Case.Fo.v_nudge[k] - gm.V.values[k]) / Case.Fo.nudge_tau
            if Case.Fo.apply_subsidence
                # Apply large-scale subsidence tendencies
                gm_H_subsidence_k = -∇H * Case.Fo.subsidence[k]
                gm_QT_subsidence_k = -∇q_tot * Case.Fo.subsidence[k]
            else
                gm_H_subsidence_k = 0.0
                gm_QT_subsidence_k = 0.0
            end

            gm.H.tendencies[k] += H_horz_adv + H_nudge + H_fluc + gm_H_subsidence_k
            gm.QT.tendencies[k] +=
                Case.Fo.dqtdt_hadv[k] + Case.Fo.dqtdt_nudge[k] + gm_QT_subsidence_k + Case.Fo.dqtdt_fluc[k]

            gm.U.tendencies[k] += gm_U_nudge_k
            gm.V.tendencies[k] += gm_V_nudge_k
        end
        gm.QT.tendencies[k] +=
            edmf.UpdThermo.qt_tendency_rain_formation_tot[k] +
            edmf.EnvThermo.qt_tendency_rain_formation[k] +
            edmf.RainPhys.qt_tendency_rain_evap[k]
        gm.H.tendencies[k] +=
            edmf.UpdThermo.θ_liq_ice_tendency_rain_formation_tot[k] +
            edmf.EnvThermo.θ_liq_ice_tendency_rain_formation[k] +
            edmf.RainPhys.θ_liq_ice_tendency_rain_evap[k]
    end

    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    edmf.massflux_h .= 0.0
    edmf.massflux_qt .= 0.0
    # Compute the mass flux and associated scalar fluxes
    @inbounds for i in 1:(up.n_updrafts)
        edmf.m[i, kf_surf] = 0.0
        a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
        @inbounds for k in real_face_indices(grid)
            a_up = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            a_en = interpc2f(ae, grid, k; a_up_bcs...)
            edmf.m[i, k] = ρ0_f[k] * a_up * a_en * (prog_up_f[i].w[k] - en.W.values[k])
        end
    end

    @inbounds for k in real_face_indices(grid)
        edmf.massflux_h[k] = 0.0
        edmf.massflux_qt[k] = 0.0
        # We know that, since W = 0 at z = 0, m = 0 also, and
        # therefore θ_liq_ice / q_tot values do not matter
        m_bcs = (; bottom = SetValue(0), top = SetValue(0))
        h_en_f = interpc2f(en.H.values, grid, k; m_bcs...)
        qt_en_f = interpc2f(en.QT.values, grid, k; m_bcs...)
        @inbounds for i in 1:(up.n_updrafts)
            h_up_f = interpc2f(prog_up[i].θ_liq_ice, grid, k; m_bcs...)
            qt_up_f = interpc2f(prog_up[i].q_tot, grid, k; m_bcs...)
            edmf.massflux_h[k] += edmf.m[i, k] * (h_up_f - h_en_f)
            edmf.massflux_qt[k] += edmf.m[i, k] * (qt_up_f - qt_en_f)
        end
    end

    # Compute the  mass flux tendencies
    # Adjust the values of the grid mean variables
    @inbounds for k in real_center_indices(grid)
        mf_tend_h_dual = dual_faces(edmf.massflux_h, grid, k)
        mf_tend_qt_dual = dual_faces(edmf.massflux_qt, grid, k)

        ∇mf_tend_h = ∇f2c(mf_tend_h_dual, grid, k)
        ∇mf_tend_qt = ∇f2c(mf_tend_qt_dual, grid, k)

        mf_tend_h = -∇mf_tend_h * α0_c[k]
        mf_tend_qt = -∇mf_tend_qt * α0_c[k]

        # Prepare the output
        edmf.massflux_tendency_h[k] = mf_tend_h
        edmf.massflux_tendency_qt[k] = mf_tend_qt
        gm.H.tendencies[k] += mf_tend_h
        gm.QT.tendencies[k] += mf_tend_qt
    end

    aeKHq_tot_bc = Case.Sur.rho_qtflux / edmf.ae[kc_surf]
    aeKHθ_liq_ice_bc = Case.Sur.rho_hflux / edmf.ae[kc_surf]
    aeKHu_bc = Case.Sur.rho_uflux / edmf.ae[kc_surf]
    aeKHv_bc = Case.Sur.rho_vflux / edmf.ae[kc_surf]
    @inbounds for k in real_center_indices(grid)
        aeKH_q_tot_cut = dual_faces(edmf.diffusive_flux_qt, grid, k)
        ∇aeKH_q_tot = ∇f2c(aeKH_q_tot_cut, grid, k; bottom = SetValue(aeKHq_tot_bc), top = SetValue(0))
        gm.QT.tendencies[k] += -α0_c[k] * ae[k] * ∇aeKH_q_tot

        aeKH_θ_liq_ice_cut = dual_faces(edmf.diffusive_flux_h, grid, k)
        ∇aeKH_θ_liq_ice = ∇f2c(aeKH_θ_liq_ice_cut, grid, k; bottom = SetValue(aeKHθ_liq_ice_bc), top = SetValue(0))
        gm.H.tendencies[k] += -α0_c[k] * ae[k] * ∇aeKH_θ_liq_ice

        aeKM_u_cut = dual_faces(edmf.diffusive_flux_u, grid, k)
        ∇aeKM_u = ∇f2c(aeKM_u_cut, grid, k; bottom = SetValue(aeKHu_bc), top = SetValue(0))
        gm.U.tendencies[k] += -α0_c[k] * ae[k] * ∇aeKM_u

        aeKM_v_cut = dual_faces(edmf.diffusive_flux_v, grid, k)
        ∇aeKM_v = ∇f2c(aeKM_v_cut, grid, k; bottom = SetValue(aeKHv_bc), top = SetValue(0))
        gm.V.tendencies[k] += -α0_c[k] * ae[k] * ∇aeKM_v
    end
end

function compute_diffusive_fluxes(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    GMV::GridMeanVariables,
    Case::CasesBase,
    TS::TimeStepping,
    param_set,
)
    ρ0_f = face_ref_state(state).ρ0
    aux_tc = center_aux_tc(state)
    edmf.ae .= 1 .- vec(aux_tc.bulk.area) # area of environment
    KM = diffusivity_m(edmf).values
    KH = diffusivity_h(edmf).values
    aeKM = edmf.ae .* KM
    aeKH = edmf.ae .* KH
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    kf_surf = kf_surface(grid)
    aeKM_bcs = (; bottom = SetValue(aeKM[kc_surf]), top = SetValue(aeKM[kc_toa]))
    aeKH_bcs = (; bottom = SetValue(aeKH[kc_surf]), top = SetValue(aeKH[kc_toa]))

    @inbounds for k in real_face_indices(grid)
        edmf.rho_ae_KH[k] = interpc2f(aeKH, grid, k; aeKH_bcs...) * ρ0_f[k]
        edmf.rho_ae_KM[k] = interpc2f(aeKM, grid, k; aeKM_bcs...) * ρ0_f[k]
    end

    aeKHq_tot_bc = -Case.Sur.rho_qtflux / edmf.ae[kc_surf] / edmf.rho_ae_KH[kc_surf]
    aeKHθ_liq_ice_bc = -Case.Sur.rho_hflux / edmf.ae[kc_surf] / edmf.rho_ae_KH[kc_surf]
    aeKHu_bc = -Case.Sur.rho_uflux / edmf.ae[kc_surf] / edmf.rho_ae_KM[kc_surf]
    aeKHv_bc = -Case.Sur.rho_vflux / edmf.ae[kc_surf] / edmf.rho_ae_KM[kc_surf]

    @inbounds for k in real_face_indices(grid)
        q_dual = dual_centers(edmf.EnvVar.QT.values, grid, k)
        ∇q_tot_f = ∇c2f(q_dual, grid, k; bottom = SetGradient(aeKHq_tot_bc), top = SetGradient(0))
        edmf.diffusive_flux_qt[k] = -edmf.rho_ae_KH[k] * ∇q_tot_f

        θ_liq_ice_dual = dual_centers(edmf.EnvVar.H.values, grid, k)
        ∇θ_liq_ice_f = ∇c2f(θ_liq_ice_dual, grid, k; bottom = SetGradient(aeKHθ_liq_ice_bc), top = SetGradient(0))
        edmf.diffusive_flux_h[k] = -edmf.rho_ae_KH[k] * ∇θ_liq_ice_f

        u_dual = dual_centers(GMV.U.values, grid, k)
        ∇u_f = ∇c2f(u_dual, grid, k; bottom = SetGradient(aeKHu_bc), top = SetGradient(0))
        edmf.diffusive_flux_u[k] = -edmf.rho_ae_KM[k] * ∇u_f

        v_dual = dual_centers(GMV.V.values, grid, k)
        ∇v_f = ∇c2f(v_dual, grid, k; bottom = SetGradient(aeKHv_bc), top = SetGradient(0))
        edmf.diffusive_flux_v[k] = -edmf.rho_ae_KM[k] * ∇v_f
    end
    return
end

# Perform the update of the scheme
function update(edmf::EDMF_PrognosticTKE, grid, state, GMV::GridMeanVariables, Case::CasesBase, TS::TimeStepping)

    gm = GMV
    up = edmf.UpdVar
    en = edmf.EnvVar
    param_set = parameter_set(gm)
    up_thermo = edmf.UpdThermo
    en_thermo = edmf.EnvThermo
    n_updrafts = up.n_updrafts

    # Update aux / pre-tendencies filters. TODO: combine these into a function that minimizes traversals
    # Some of these methods should probably live in `compute_tendencies`, when written, but we'll
    # treat them as auxiliary variables for now, until we disentangle the tendency computations.
    set_updraft_surface_bc(edmf, grid, state, gm, Case)
    update_aux!(edmf, gm, grid, state, Case, param_set, TS)

    compute_rain_formation_tendencies(up_thermo, grid, state, edmf.UpdVar, edmf.Rain, TS.dt, param_set) # causes division error in dry bubble first time step
    microphysics(en_thermo, grid, state, edmf.EnvVar, edmf.Rain, TS.dt, param_set) # saturation adjustment + rain creation
    if edmf.Rain.rain_model == "clima_1m"
        compute_rain_evap_tendencies(edmf.RainPhys, grid, state, gm, TS, edmf.Rain.QR)
        compute_rain_advection_tendencies(edmf.RainPhys, grid, state, gm, TS, edmf.Rain.QR)
    end

    # compute tendencies
    compute_gm_tendencies!(edmf, grid, state, Case, gm, TS)
    compute_updraft_tendencies(edmf, grid, state, gm, TS)
    # ----------- TODO: move to compute_tendencies
    implicit_eqs = edmf.implicit_eqs
    # Matrix is the same for all variables that use the same eddy diffusivity, we can construct once and reuse
    prog_up = center_prog_updrafts(state)

    KM = diffusivity_m(edmf).values
    KH = diffusivity_h(edmf).values
    common_args = (
        grid,
        param_set,
        state,
        TS,
        KM,
        KH,
        en.W.values,
        en.TKE.values,
        up.n_updrafts,
        edmf.minimum_area,
        edmf.pressure_plume_spacing,
        edmf.frac_turb_entr,
        edmf.entr_sc,
        edmf.mixing_length,
    )

    implicit_eqs.A_TKE .= construct_tridiag_diffusion_en(common_args..., true)
    implicit_eqs.A_Hvar .= construct_tridiag_diffusion_en(common_args..., false)
    implicit_eqs.A_QTvar .= construct_tridiag_diffusion_en(common_args..., false)
    implicit_eqs.A_HQTcov .= construct_tridiag_diffusion_en(common_args..., false)

    implicit_eqs.b_TKE .= en_diffusion_tendencies(grid, state, TS, en.TKE, n_updrafts)
    implicit_eqs.b_Hvar .= en_diffusion_tendencies(grid, state, TS, en.Hvar, n_updrafts)
    implicit_eqs.b_QTvar .= en_diffusion_tendencies(grid, state, TS, en.QTvar, n_updrafts)
    implicit_eqs.b_HQTcov .= en_diffusion_tendencies(grid, state, TS, en.HQTcov, n_updrafts)
    # -----------

    ###
    ### update
    ###
    update_updraft(edmf, grid, state, gm, TS)
    if edmf.Rain.rain_model == "clima_1m"
        update_rain(edmf.Rain, grid, state, up_thermo, en_thermo, edmf.RainPhys, TS)
    end

    en.TKE.values .= tridiag_solve(implicit_eqs.b_TKE, implicit_eqs.A_TKE)
    en.Hvar.values .= tridiag_solve(implicit_eqs.b_Hvar, implicit_eqs.A_Hvar)
    en.QTvar.values .= tridiag_solve(implicit_eqs.b_QTvar, implicit_eqs.A_QTvar)
    en.HQTcov.values .= tridiag_solve(implicit_eqs.b_HQTcov, implicit_eqs.A_HQTcov)
    @inbounds for k in real_center_indices(grid)
        gm.U.values[k] += gm.U.tendencies[k] * TS.dt
        gm.V.values[k] += gm.V.tendencies[k] * TS.dt
        gm.H.values[k] += gm.H.tendencies[k] * TS.dt
        gm.QT.values[k] += gm.QT.tendencies[k] * TS.dt
    end

    ###
    ### set values
    ###
    @inbounds for k in real_center_indices(grid)
        en.TKE.values[k] = max(en.TKE.values[k], 0.0)
        en.Hvar.values[k] = max(en.Hvar.values[k], 0.0)
        en.QTvar.values[k] = max(en.QTvar.values[k], 0.0)
        en.HQTcov.values[k] = max(en.HQTcov.values[k], -sqrt(en.Hvar.values[k] * en.QTvar.values[k]))
        en.HQTcov.values[k] = min(en.HQTcov.values[k], sqrt(en.Hvar.values[k] * en.QTvar.values[k]))
    end

    # set values
    set_values_with_new(edmf.UpdVar, grid, state)
    return
end

function set_updraft_surface_bc(edmf::EDMF_PrognosticTKE, grid, state, GMV::GridMeanVariables, Case::CasesBase)
    kc_surf = kc_surface(grid)

    Δzi = grid.Δzi
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]
    qt_var = get_surface_variance(Case.Sur.rho_qtflux * α0LL, Case.Sur.rho_qtflux * α0LL, ustar, zLL, oblength)
    h_var = get_surface_variance(Case.Sur.rho_hflux * α0LL, Case.Sur.rho_hflux * α0LL, ustar, zLL, oblength)

    if Case.Sur.bflux > 0.0
        a_total = edmf.surface_area
        edmf.entr_surface_bc = 2.0 * Δzi
        edmf.detr_surface_bc = 0.0
    else
        # a_total = edmf.surface_area
        a_total = edmf.minimum_area * 0.9
        edmf.entr_surface_bc = 0.0
        edmf.detr_surface_bc = 2.0 * Δzi
    end

    a_ = a_total / edmf.n_updrafts
    @inbounds for i in 1:(edmf.n_updrafts)
        surface_scalar_coeff = percentile_bounds_mean_norm(1.0 - a_total + i * a_, 1.0 - a_total + (i + 1) * a_, 1000)
        edmf.area_surface_bc[i] = a_
        edmf.w_surface_bc[i] = 0.0
        edmf.h_surface_bc[i] = (GMV.H.values[kc_surf] + surface_scalar_coeff * sqrt(h_var))
        edmf.qt_surface_bc[i] = (GMV.QT.values[kc_surf] + surface_scalar_coeff * sqrt(qt_var))
    end
    return
end

function reset_surface_covariance(edmf::EDMF_PrognosticTKE, grid, state, GMV, Case::CasesBase)
    flux1 = Case.Sur.rho_hflux
    flux2 = Case.Sur.rho_qtflux
    kc_surf = kc_surface(grid)
    zLL = grid.zc[kc_surf]
    ustar = Case.Sur.ustar
    oblength = Case.Sur.obukhov_length
    α0LL = center_ref_state(state).α0[kc_surf]

    gm = GMV
    up = edmf.UpdVar
    en = edmf.EnvVar

    en.TKE.values[kc_surf] = get_surface_tke(Case.Sur.ustar, grid.zc[kc_surf], Case.Sur.obukhov_length)
    get_GMV_CoVar(edmf, grid, state, en.W, en.W, en.TKE, gm.W.values, gm.W.values, gm.TKE.values, :w)

    en.Hvar.values[kc_surf] = get_surface_variance(flux1 * α0LL, flux1 * α0LL, ustar, zLL, oblength)
    en.QTvar.values[kc_surf] = get_surface_variance(flux2 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    en.HQTcov.values[kc_surf] = get_surface_variance(flux1 * α0LL, flux2 * α0LL, ustar, zLL, oblength)
    return
end

# Note: this assumes all variables are defined on half levels not full levels (i.e. phi, psi are not w)
# if covar_e.name is not "tke".
function get_GMV_CoVar(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    ϕ_en::EnvironmentVariable,
    ψ_en::EnvironmentVariable,
    covar_e::EnvironmentVariable_2m,
    ϕ_gm,
    ψ_gm,
    gmv_covar,
    ϕ_up_sym::Symbol,
    ψ_up_sym::Symbol = ϕ_up_sym,
)

    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    is_tke = covar_e.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)

    if is_tke
        @inbounds for k in real_center_indices(grid)
            ϕ_en_dual = dual_faces(ϕ_en.values, grid, k)
            ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
            ψ_en_dual = dual_faces(ψ_en.values, grid, k)
            ψ_gm_dual = dual_faces(ψ_gm, grid, k)
            Δϕ_dual = ϕ_en_dual .- ϕ_gm_dual
            Δψ_dual = ψ_en_dual .- ψ_gm_dual
            Δϕ = interpf2c(Δϕ_dual, grid, k)
            Δψ = interpf2c(Δψ_dual, grid, k)

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e.values[k]
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up_var = getproperty(prog_up_f[i], ϕ_up_sym)
                ψ_up_var = getproperty(prog_up_f[i], ψ_up_sym)
                ϕ_up_dual = dual_faces(ϕ_up_var, grid, k)
                ϕ_gm_dual = dual_faces(ϕ_gm, grid, k)
                ψ_up_dual = dual_faces(ψ_up_var, grid, k)
                ψ_gm_dual = dual_faces(ψ_gm, grid, k)
                Δϕ_dual = ϕ_up_dual .- ϕ_gm_dual
                Δψ_dual = ψ_up_dual .- ψ_gm_dual
                Δϕ = interpf2c(Δϕ_dual, grid, k)
                Δψ = interpf2c(Δψ_dual, grid, k)
                gmv_covar[k] += tke_factor * prog_up[i].area[k] * Δϕ * Δψ
            end
        end
    else

        @inbounds for k in real_center_indices(grid)
            Δϕ = ϕ_en.values[k] - ϕ_gm[k]
            Δψ = ψ_en.values[k] - ψ_gm[k]

            gmv_covar[k] = tke_factor * ae[k] * Δϕ * Δψ + ae[k] * covar_e.values[k]
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up_var = getproperty(prog_up[i], ϕ_up_sym)
                ψ_up_var = getproperty(prog_up[i], ψ_up_sym)
                Δϕ = ϕ_up_var[k] - ϕ_gm[k]
                Δψ = ψ_up_var[k] - ψ_gm[k]
                gmv_covar[k] += tke_factor * prog_up[i].area[k] * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_pressure_plume_spacing(edmf::EDMF_PrognosticTKE, param_set)

    H_up_min = CPEDMF.H_up_min(param_set)
    @inbounds for i in 1:(edmf.n_updrafts)
        edmf.pressure_plume_spacing[i] =
            max(edmf.aspect_ratio * edmf.UpdVar.updraft_top[i], H_up_min * edmf.aspect_ratio)
    end
    return
end

function compute_updraft_tendencies(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, TS::TimeStepping)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dti_ = 1.0 / TS.dt
    Δt = TS.dt

    up = edmf.UpdVar
    up_thermo = edmf.UpdThermo
    en = edmf.EnvVar
    prog_up = center_prog_updrafts(state)
    aux_up = center_aux_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    ρ_0_c = center_ref_state(state).ρ0
    ρ_0_f = face_ref_state(state).ρ0
    au_lim = edmf.max_area

    @inbounds for i in 1:(up.n_updrafts)
        @inbounds for k in real_center_indices(grid)
            tendencies_up[i].area[k] = 0
            tendencies_up[i].θ_liq_ice[k] = 0
            tendencies_up[i].q_tot[k] = 0
        end
        @inbounds for k in real_face_indices(grid)
            tendencies_up_f[i].w[k] = 0
        end
    end

    @inbounds for i in 1:(up.n_updrafts)
        edmf.entr_sc[i, kc_surf] = edmf.entr_surface_bc
        edmf.detr_sc[i, kc_surf] = edmf.detr_surface_bc
    end

    # Solve for updraft area fraction
    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            adv = upwind_advection_area(ρ_0_c, prog_up[i].area, prog_up_f[i].w, grid, k)

            a_up_c = prog_up[i].area[k]
            entr_term = a_up_c * w_up_c * (edmf.entr_sc[i, k])
            detr_term = a_up_c * w_up_c * (-edmf.detr_sc[i, k])
            tendencies_up[i].area[k] = adv + entr_term + detr_term
            a_up_candidate = a_up_c + Δt * tendencies_up[i].area[k]
            if a_up_candidate > au_lim
                a_up_div = a_up_c > 0.0 ? a_up_c : au_lim
                a_up_candidate = au_lim
                edmf.detr_sc[i, k] = (((au_lim - a_up_c) * dti_ - adv - entr_term) / (-a_up_div * w_up_c))
                tendencies_up[i].area[k] = (a_up_candidate - a_up_c) * dti_
            end
        end
    end

    entr_w_c = edmf.entr_sc .+ edmf.frac_turb_entr
    detr_w_c = edmf.detr_sc .+ edmf.frac_turb_entr

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            m_k = (ρ_0_c[k] * prog_up[i].area[k] * w_up_c)

            adv = upwind_advection_scalar(ρ_0_c, prog_up[i].area, prog_up_f[i].w, prog_up[i].θ_liq_ice, grid, k)
            entr = entr_w_c[i, k] * en.H.values[k]
            detr = detr_w_c[i, k] * prog_up[i].θ_liq_ice[k]
            rain = ρ_0_c[k] * up_thermo.θ_liq_ice_tendency_rain_formation[i, k]
            tendencies_up[i].θ_liq_ice[k] = -adv + m_k * (entr - detr) + rain

            adv = upwind_advection_scalar(ρ_0_c, prog_up[i].area, prog_up_f[i].w, prog_up[i].q_tot, grid, k)
            entr = entr_w_c[i, k] * en.QT.values[k]
            detr = detr_w_c[i, k] * prog_up[i].q_tot[k]
            rain = ρ_0_c[k] * up_thermo.qt_tendency_rain_formation[i, k]
            tendencies_up[i].q_tot[k] = -adv + m_k * (entr - detr) + rain
        end
    end

    # Solve for updraft velocity
    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            a_k = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
            # We know that, since W = 0 at z = 0, these BCs should
            # not matter in the end:
            entr_w = interpc2f(entr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            detr_w = interpc2f(detr_w_c, grid, k, i; bottom = SetValue(0), top = SetValue(0))
            B_k = interpc2f(aux_up[i].buoy, grid, k; bottom = SetValue(0), top = SetValue(0))

            adv = upwind_advection_velocity(ρ_0_f, prog_up[i].area, prog_up_f[i].w, grid, k; a_up_bcs)
            exch = (ρ_0_f[k] * a_k * prog_up_f[i].w[k] * (entr_w * en.W.values[k] - detr_w * prog_up_f[i].w[k]))
            buoy = ρ_0_f[k] * a_k * B_k
            tendencies_up_f[i].w[k] = -adv + exch + buoy + edmf.nh_pressure[i, k]
        end
    end
    return
end


function update_updraft(edmf::EDMF_PrognosticTKE, grid, state, gm::GridMeanVariables, TS::TimeStepping)
    param_set = parameter_set(gm)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    dti_ = 1.0 / TS.dt
    Δt = TS.dt

    up = edmf.UpdVar
    up_thermo = edmf.UpdThermo
    en = edmf.EnvVar
    tendencies_up = center_tendencies_updrafts(state)
    tendencies_up_f = face_tendencies_updrafts(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    a_up_new = up.Area.new
    w_up_new = up.W.new
    ρ_0_c = center_ref_state(state).ρ0
    ρ_0_f = face_ref_state(state).ρ0

    @inbounds for i in 1:(up.n_updrafts)
        w_up_new[i, kf_surf] = edmf.w_surface_bc[i]
        a_up_new[i, kc_surf] = edmf.area_surface_bc[i]
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            is_surface_center(grid, k) && continue
            a_up_c = prog_up[i].area[k]
            a_up_candidate = max(a_up_c + Δt * tendencies_up[i].area[k], 0)
            a_up_new[i, k] = a_up_candidate
        end
    end

    @inbounds for k in real_face_indices(grid)
        is_surface_face(grid, k) && continue
        @inbounds for i in 1:(up.n_updrafts)
            a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
            anew_k = interpc2f(a_up_new, grid, k, i; a_up_bcs...)
            if anew_k >= edmf.minimum_area
                a_k = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
                w_up_new[i, k] =
                    (ρ_0_f[k] * a_k * prog_up_f[i].w[k] + Δt * tendencies_up_f[i].w[k]) / (ρ_0_f[k] * anew_k)

                w_up_new[i, k] = max(w_up_new[i, k], 0)
                # TODO: remove a_up_new from this loop.
                if w_up_new[i, k] <= 0.0
                    if !(k.i > size(a_up_new, 2))
                        a_up_new[i, k] = 0
                    end
                end
            else
                w_up_new[i, k] = 0
                if !(k.i > size(a_up_new, 2))
                    a_up_new[i, k] = 0
                end
            end
        end
    end

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            if is_surface_center(grid, k)
                # at the surface
                if a_up_new[i, k] >= edmf.minimum_area
                    up.H.new[i, k] = edmf.h_surface_bc[i]
                    up.QT.new[i, k] = edmf.qt_surface_bc[i]
                else
                    up.H.new[i, k] = gm.H.values[k]
                    up.QT.new[i, k] = gm.QT.values[k]
                end
                continue
            end

            if a_up_new[i, k] >= edmf.minimum_area
                up.H.new[i, k] =
                    (ρ_0_c[k] * prog_up[i].area[k] * prog_up[i].θ_liq_ice[k] + Δt * tendencies_up[i].θ_liq_ice[k]) /
                    (ρ_0_c[k] * a_up_new[i, k])

                up.QT.new[i, k] = max(
                    (ρ_0_c[k] * prog_up[i].area[k] * prog_up[i].q_tot[k] + Δt * tendencies_up[i].q_tot[k]) /
                    (ρ_0_c[k] * a_up_new[i, k]),
                    0.0,
                )

            else
                up.H.new[i, k] = gm.H.values[k]
                up.QT.new[i, k] = gm.QT.values[k]
            end
        end
    end

    return
end

function initialize_covariance(edmf::EDMF_PrognosticTKE, grid, state, gm, Case::CasesBase)

    kc_surf = kc_surface(grid)
    en = edmf.EnvVar

    en.TKE.values .= gm.TKE.values

    reset_surface_covariance(edmf, grid, state, gm, Case)
    gm.Hvar.values .= gm.Hvar.values[kc_surf] .* gm.TKE.values
    gm.QTvar.values .= gm.QTvar.values[kc_surf] .* gm.TKE.values
    gm.HQTcov.values .= gm.HQTcov.values[kc_surf] .* gm.TKE.values

    en.Hvar.values .= gm.Hvar.values
    en.QTvar.values .= gm.QTvar.values
    en.HQTcov.values .= gm.HQTcov.values
    return
end

function compute_covariance_shear(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    gm::GridMeanVariables,
    Covar::EnvironmentVariable_2m,
    EnvVar1,
    EnvVar2,
)

    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area # area of environment
    KH = diffusivity_h(edmf).values
    ρ0_c = center_ref_state(state).ρ0
    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    k_eddy = is_tke ? diffusivity_m(edmf).values : diffusivity_h(edmf).values

    if is_tke
        @inbounds for k in real_center_indices(grid)
            v_cut = ccut(gm.V.values, grid, k)
            ∇v = c∇(v_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            u_cut = ccut(gm.U.values, grid, k)
            ∇u = c∇(u_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_dual = dual_faces(EnvVar2, grid, k)
            var1_dual = dual_faces(EnvVar1, grid, k)

            ∇var2 = ∇f2c(var2_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))
            ∇var1 = ∇f2c(var1_dual, grid, k; bottom = SetValue(0), top = SetGradient(0))

            Covar.shear[k] = tke_factor * 2 * (ρ0_c[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2 + ∇u^2 + ∇v^2))
        end
    else
        @inbounds for k in real_center_indices(grid)
            # Defined correctly only for covariance between half-level variables.
            var1_cut = ccut(EnvVar1, grid, k)
            ∇var1 = c∇(var1_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            var2_cut = ccut(EnvVar2, grid, k)
            ∇var2 = c∇(var2_cut, grid, k; bottom = Extrapolate(), top = SetGradient(0))

            Covar.shear[k] = tke_factor * 2 * (ρ0_c[k] * ae[k] * k_eddy[k] * (∇var1 * ∇var2))
        end
    end
    return
end

function compute_covariance_interdomain_src(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    ϕ_en::EnvironmentVariable,
    ψ_en::EnvironmentVariable,
    Covar::EnvironmentVariable_2m,
    ϕ_var::Symbol,
    ψ_var::Symbol = ϕ_var,
)

    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    prog_up_c = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    prog_up = is_tke ? prog_up_f : prog_up_c
    if is_tke
        @inbounds for k in real_center_indices(grid)
            Covar.interdomain[k] = 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up = getproperty(prog_up[i], ϕ_var)
                ψ_up = getproperty(prog_up[i], ψ_var)
                Δϕ = interpf2c(ϕ_up, grid, k) - interpf2c(ϕ_en.values, grid, k)
                Δψ = interpf2c(ψ_up, grid, k) - interpf2c(ψ_en.values, grid, k)

                Covar.interdomain[k] += tke_factor * prog_up_c[i].area[k] * (1.0 - prog_up_c[i].area[k]) * Δϕ * Δψ
            end
        end
    else
        @inbounds for k in real_center_indices(grid)
            Covar.interdomain[k] = 0.0
            @inbounds for i in 1:(edmf.n_updrafts)
                ϕ_up = getproperty(prog_up[i], ϕ_var)
                ψ_up = getproperty(prog_up[i], ψ_var)
                Δϕ = ϕ_up[k] - ϕ_en.values[k]
                Δψ = ψ_up[k] - ψ_en.values[k]
                Covar.interdomain[k] += tke_factor * prog_up_c[i].area[k] * (1.0 - prog_up_c[i].area[k]) * Δϕ * Δψ
            end
        end
    end
    return
end

function compute_covariance_entr(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    Covar::EnvironmentVariable_2m,
    EnvVar1::EnvironmentVariable,
    EnvVar2::EnvironmentVariable,
    GmvVar1::VariablePrognostic,
    GmvVar2::VariablePrognostic,
    var1::Symbol,
    var2::Symbol = var1,
)

    ρ_0_c = center_ref_state(state).ρ0

    is_tke = Covar.name == "tke"
    tke_factor = is_tke ? 0.5 : 1
    prog_up_c = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    prog_up = is_tke ? prog_up_f : prog_up_c

    @inbounds for k in real_center_indices(grid)
        Covar.entr_gain[k] = 0.0
        Covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(edmf.n_updrafts)
            a_up = prog_up_c[i].area[k]
            if a_up > edmf.minimum_area
                R_up = edmf.pressure_plume_spacing[i]
                up_var_1 = getproperty(prog_up[i], var1)
                up_var_2 = getproperty(prog_up[i], var2)
                updvar1 = is_tke ? interpf2c(up_var_1, grid, k) : up_var_1[k]
                updvar2 = is_tke ? interpf2c(up_var_2, grid, k) : up_var_2[k]
                envvar1 = is_tke ? interpf2c(EnvVar1.values, grid, k) : EnvVar1.values[k]
                envvar2 = is_tke ? interpf2c(EnvVar2.values, grid, k) : EnvVar2.values[k]
                gmvvar1 = is_tke ? interpf2c(GmvVar1.values, grid, k) : GmvVar1.values[k]
                gmvvar2 = is_tke ? interpf2c(GmvVar2.values, grid, k) : GmvVar2.values[k]

                eps_turb = edmf.frac_turb_entr[i, k]

                w_u = interpf2c(prog_up_f[i].w, grid, k)
                dynamic_entr =
                    tke_factor *
                    ρ_0_c[k] *
                    a_up *
                    abs(w_u) *
                    edmf.detr_sc[i, k] *
                    (updvar1 - envvar1) *
                    (updvar2 - envvar2)
                turbulent_entr =
                    tke_factor *
                    ρ_0_c[k] *
                    a_up *
                    abs(w_u) *
                    eps_turb *
                    ((envvar1 - gmvvar1) * (updvar2 - envvar2) + (envvar2 - gmvvar2) * (updvar1 - envvar1))
                Covar.entr_gain[k] += dynamic_entr + turbulent_entr
                Covar.detr_loss[k] +=
                    tke_factor * ρ_0_c[k] * a_up * abs(w_u) * (edmf.entr_sc[i, k] + eps_turb) * Covar.values[k]
            end
        end
    end

    return
end

function compute_covariance_detr(edmf::EDMF_PrognosticTKE, grid, state, Covar::EnvironmentVariable_2m)
    up = edmf.UpdVar
    ρ0_c = center_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    @inbounds for k in real_center_indices(grid)
        Covar.detr_loss[k] = 0.0
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            Covar.detr_loss[k] += prog_up[i].area[k] * abs(w_up_c) * edmf.entr_sc[i, k]
        end
        Covar.detr_loss[k] *= ρ0_c[k] * Covar.values[k]
    end
    return
end

function compute_covariance_dissipation(edmf::EDMF_PrognosticTKE, grid, state, Covar::EnvironmentVariable_2m, param_set)
    en = edmf.EnvVar
    c_d = CPEDMF.c_d(param_set)
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    ρ0_c = center_ref_state(state).ρ0

    @inbounds for k in real_center_indices(grid)
        Covar.dissipation[k] = (
            ρ0_c[k] * ae[k] * Covar.values[k] * max(en.TKE.values[k], 0)^0.5 / max(edmf.mixing_length[k], 1.0e-3) * c_d
        )
    end
    return
end

function en_diffusion_tendencies(grid::Grid, state, TS, covar, n_updrafts)
    dti = TS.dti
    b = center_field(grid)
    ρ0_c = center_ref_state(state).ρ0
    prog_up = center_prog_updrafts(state)

    ae = center_field(grid)

    @inbounds for k in real_center_indices(grid)
        ae[k] = 1 .- sum(ntuple(i -> prog_up[i].area[k], n_updrafts))
    end

    kc_surf = kc_surface(grid)
    covar_surf = covar.values[kc_surf]

    @inbounds for k in real_center_indices(grid)
        if is_surface_center(grid, k)
            b[k] = covar_surf
        else
            b[k] = (
                ρ0_c[k] * ae[k] * covar.values[k] * dti +
                covar.press[k] +
                covar.buoy[k] +
                covar.shear[k] +
                covar.entr_gain[k] +
                covar.rain_src[k]
            )
        end
    end

    return b
end

function GMV_third_m(
    edmf::EDMF_PrognosticTKE,
    grid,
    state,
    Gmv_third_m::VariableDiagnostic,
    env_covar::EnvironmentVariable_2m,
    env_mean::EnvironmentVariable,
    up_var::Symbol,
)

    up = edmf.UpdVar
    en = edmf.EnvVar
    aux_tc = center_aux_tc(state)
    ae = 1 .- aux_tc.bulk.area
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    is_tke = env_covar.name == "tke"

    @inbounds for k in real_center_indices(grid)
        mean_en = is_tke ? interpf2c(env_mean.values, grid, k) : env_mean.values[k]
        GMVv_ = ae[k] * mean_en
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(prog_up_f[i], up_var) : getproperty(prog_up[i], up_var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVv_ += prog_up[i].area[k] * mean_up
        end

        # TODO: report bug: i used outside of scope.
        # This is only valid (assuming correct) for 1
        # updraft.
        i_last = last(1:(up.n_updrafts))
        if is_tke
            w_bcs = (; bottom = SetValue(0), top = SetValue(0))
            w_en_dual = dual_faces(en.W.values, grid, k)
            ∇w_en = ∇f2c(w_en_dual, grid, k; w_bcs...)
            Envcov_ = -edmf.horiz_K_eddy[i_last, k] * ∇w_en
        else
            Envcov_ = env_covar.values[k]
        end

        Upd_cubed = 0.0
        GMVcov_ = ae[k] * (Envcov_ + (mean_en - GMVv_)^2)
        @inbounds for i in 1:(up.n_updrafts)
            var_up = is_tke ? getproperty(prog_up_f[i], up_var) : getproperty(prog_up[i], up_var)
            mean_up = is_tke ? interpf2c(var_up, grid, k) : var_up[k]
            GMVcov_ += prog_up[i].area[k] * (mean_up - GMVv_)^2
            Upd_cubed += prog_up[i].area[k] * mean_up^3
        end

        if is_surface_center(grid, k)
            Gmv_third_m.values[k] = 0.0 # this is here as first value is biased with BC area fraction
        else
            Gmv_third_m.values[k] =
                Upd_cubed + ae[k] * (mean_en^3 + 3 * mean_en * Envcov_) - GMVv_^3 - 3 * GMVcov_ * GMVv_
        end
    end
    return
end

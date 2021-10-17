function update_aux!(edmf, gm, grid, state, Case, param_set, TS)
    #####
    ##### Unpack common variables
    #####
    kc_surf = kc_surface(grid)
    kc_toa = kc_top_of_atmos(grid)
    up = edmf.UpdVar
    en = edmf.EnvVar
    en_thermo = edmf.EnvThermo
    ρ0_f = face_ref_state(state).ρ0
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    α0_c = center_ref_state(state).α0
    g = CPP.grav(param_set)
    c_m = CPEDMF.c_m(param_set)
    KM = center_aux_tc(state).KM
    KH = center_aux_tc(state).KH
    surface = Case.Sur
    obukhov_length = surface.obukhov_length
    FT = eltype(grid)
    prog_up = center_prog_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_gm = center_aux_grid_mean(state)
    prog_up_f = face_prog_updrafts(state)
    aux_up_f = face_aux_tc(state)
    aux_tc = center_aux_tc(state)
    prog_en = center_prog_environment(state)
    aux_en_2m = center_aux_environment_2m(state)

    for k in real_center_indices(grid)
        aux_tc.bulk.area[k] = sum(ntuple(i -> prog_up[i].area[k], up.n_updrafts))
    end

    #####
    ##### diagnose_GMV_moments
    #####
    #! format: off
    get_GMV_CoVar(edmf, grid, state, en.H, en.H, :Hvar, :θ_liq_ice)
    get_GMV_CoVar(edmf, grid, state, en.QT, en.QT, :QTvar, :q_tot)
    get_GMV_CoVar(edmf, grid, state, en.H, en.QT, :HQTcov, :θ_liq_ice, :q_tot)
    GMV_third_m(edmf, grid, state, :Hvar, en.H, :θ_liq_ice, :H_third_m)
    GMV_third_m(edmf, grid, state, :QTvar, en.QT, :q_tot, :QT_third_m)
    GMV_third_m(edmf, grid, state, :tke, en.W, :w, :W_third_m)
    #! format: on

    #####
    ##### decompose_environment
    #####
    # Find values of environmental variables by subtracting updraft values from grid mean values
    # whichvals used to check which substep we are on--correspondingly use "gm.SomeVar" (last timestep value)
    # first make sure the "bulkvalues" of the updraft variables are updated
    @inbounds for k in real_face_indices(grid)
        aux_up_f.bulk.w[k] = 0
        a_bulk_bcs = (; bottom = SetValue(sum(edmf.area_surface_bc)), top = SetZeroGradient())
        a_bulk_f = interpc2f(aux_tc.bulk.area, grid, k; a_bulk_bcs...)
        if a_bulk_f > 1.0e-20
            @inbounds for i in 1:(up.n_updrafts)
                a_up_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetZeroGradient())
                a_up_f = interpc2f(prog_up[i].area, grid, k; a_up_bcs...)
                aux_up_f.bulk.w[k] += a_up_f * prog_up_f[i].w[k] / a_bulk_f
            end
        end
        # Assuming gm.W = 0!
        en.W.values[k] = -a_bulk_f / (1 - a_bulk_f) * aux_up_f.bulk.w[k]
    end

    @inbounds for k in real_center_indices(grid)
        a_bulk_c = aux_tc.bulk.area[k]
        aux_tc.bulk.q_tot[k] = 0
        aux_tc.bulk.q_liq[k] = 0
        aux_tc.bulk.θ_liq_ice[k] = 0
        aux_tc.bulk.T[k] = 0
        aux_tc.bulk.RH[k] = 0
        aux_tc.bulk.buoy[k] = 0
        if a_bulk_c > 1.0e-20
            @inbounds for i in 1:(up.n_updrafts)
                aux_tc.bulk.q_tot[k] += prog_up[i].area[k] * prog_up[i].q_tot[k] / a_bulk_c
                aux_tc.bulk.q_liq[k] += prog_up[i].area[k] * aux_up[i].q_liq[k] / a_bulk_c
                aux_tc.bulk.θ_liq_ice[k] += prog_up[i].area[k] * prog_up[i].θ_liq_ice[k] / a_bulk_c
                aux_tc.bulk.T[k] += prog_up[i].area[k] * aux_up[i].T[k] / a_bulk_c
                aux_tc.bulk.RH[k] += prog_up[i].area[k] * aux_up[i].RH[k] / a_bulk_c
                aux_tc.bulk.buoy[k] += prog_up[i].area[k] * aux_up[i].buoy[k] / a_bulk_c
            end
        else
            aux_tc.bulk.q_tot[k] = prog_gm.q_tot[k]
            aux_tc.bulk.θ_liq_ice[k] = prog_gm.θ_liq_ice[k]
            aux_tc.bulk.RH[k] = aux_gm.RH[k]
            aux_tc.bulk.T[k] = aux_gm.T[k]
        end
        if aux_tc.bulk.q_liq[k] > 1e-8 && a_bulk_c > 1e-3
            up.cloud_fraction[k] = 1.0
        else
            up.cloud_fraction[k] = 0.0
        end

        val1 = 1 / (1 - a_bulk_c)
        val2 = a_bulk_c * val1

        en.Area.values[k] = 1 - a_bulk_c
        en.QT.values[k] = max(val1 * prog_gm.q_tot[k] - val2 * aux_tc.bulk.q_tot[k], 0) #Yair - this is here to prevent negative QT
        en.H.values[k] = val1 * prog_gm.θ_liq_ice[k] - val2 * aux_tc.bulk.θ_liq_ice[k]

        #####
        ##### saturation_adjustment
        #####

        ts_en = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], en.QT.values[k])

        en.T.values[k] = TD.air_temperature(ts_en)
        en.QL.values[k] = TD.liquid_specific_humidity(ts_en)
        rho = TD.air_density(ts_en)
        en.B.values[k] = buoyancy_c(param_set, ρ0_c[k], rho)

        update_cloud_dry(en_thermo, k, en, ts_en)
        en.RH.values[k] = TD.relative_humidity(ts_en)

        #####
        ##### buoyancy
        #####

        @inbounds for i in 1:(up.n_updrafts)
            if prog_up[i].area[k] > 0.0
                ts_up = TD.PhaseEquil_pθq(param_set, p0_c[k], prog_up[i].θ_liq_ice[k], prog_up[i].q_tot[k])
                aux_up[i].q_liq[k] = TD.liquid_specific_humidity(ts_up)
                aux_up[i].T[k] = TD.air_temperature(ts_up)
                ρ = TD.air_density(ts_up)
                aux_up[i].buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                aux_up[i].RH[k] = TD.relative_humidity(ts_up)
            elseif k > kc_surf
                if prog_up[i].area[k - 1] > 0.0 && edmf.extrapolate_buoyancy
                    qt = prog_up[i].q_tot[k - 1]
                    h = prog_up[i].θ_liq_ice[k - 1]
                    ts_up = TD.PhaseEquil_pθq(param_set, p0_c[k], h, qt)
                    ρ = TD.air_density(ts_up)
                    aux_up[i].buoy[k] = buoyancy_c(param_set, ρ0_c[k], ρ)
                    aux_up[i].RH[k] = TD.relative_humidity(ts_up)
                else
                    aux_up[i].buoy[k] = en.B.values[k]
                    aux_up[i].RH[k] = en.RH.values[k]
                end
            else
                aux_up[i].buoy[k] = en.B.values[k]
                aux_up[i].RH[k] = en.RH.values[k]
            end
        end

        aux_gm.buoy[k] = (1.0 - aux_tc.bulk.area[k]) * en.B.values[k]
        @inbounds for i in 1:(up.n_updrafts)
            aux_gm.buoy[k] += prog_up[i].area[k] * aux_up[i].buoy[k]
        end
        @inbounds for i in 1:(up.n_updrafts)
            aux_up[i].buoy[k] -= aux_gm.buoy[k]
        end
        en.B.values[k] -= aux_gm.buoy[k]
    end
    # TODO - use this inversion in free_convection_windspeed and not compute zi twice
    θ_ρ = center_field(grid)
    @inbounds for k in real_center_indices(grid)
        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], prog_gm.θ_liq_ice[k], prog_gm.q_tot[k])
        θ_ρ[k] = TD.virtual_pottemp(ts)
    end
    edmf.zi = get_inversion(param_set, θ_ρ, prog_gm.u, prog_gm.v, grid, surface.Ri_bulk_crit)

    update_surface(Case, grid, state, gm, TS, param_set)
    update_forcing(Case, grid, state, gm, TS, param_set)
    update_radiation(Case, grid, state, gm, TS, param_set)

    #####
    ##### update_GMV_diagnostics
    #####
    a_up_bulk = aux_tc.bulk.area
    @inbounds for k in real_center_indices(grid)
        aux_gm.q_liq[k] = (a_up_bulk[k] * aux_tc.bulk.q_liq[k] + (1 - a_up_bulk[k]) * en.QL.values[k])
        aux_gm.T[k] = (a_up_bulk[k] * aux_tc.bulk.T[k] + (1 - a_up_bulk[k]) * en.T.values[k])
        aux_gm.buoy[k] = (a_up_bulk[k] * aux_tc.bulk.buoy[k] + (1 - a_up_bulk[k]) * en.B.values[k])
    end
    compute_pressure_plume_spacing(edmf, param_set)

    #####
    ##### compute_updraft_closures
    #####
    upd_cloud_diagnostics(up, grid, state) # TODO: should this be moved to compute_diagnostics! ?

    @inbounds for k in real_center_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)
            # entrainment
            if prog_up[i].area[k] > 0.0
                # compute ∇m at cell centers
                a_up_c = prog_up[i].area[k]
                w_up_c = interpf2c(prog_up_f[i].w, grid, k)
                w_gm_c = interpf2c(prog_gm_f.w, grid, k)
                m = a_up_c * (w_up_c - w_gm_c)
                a_up_cut = ccut_upwind(prog_up[i].area, grid, k)
                w_up_cut = daul_f2c_upwind(prog_up_f[i].w, grid, k)
                w_gm_cut = daul_f2c_upwind(prog_gm_f.w, grid, k)
                m_cut = a_up_cut .* (w_up_cut .- w_gm_cut)
                ∇m = FT(c∇_upwind(m_cut, grid, k; bottom = SetValue(0), top = FreeBoundary()))

                w_min = 0.001

                εδ_model = MoistureDeficitEntr(;
                    q_liq_up = aux_up[i].q_liq[k],
                    q_liq_en = en.QL.values[k],
                    w_up = interpf2c(prog_up_f[i].w, grid, k),
                    w_en = interpf2c(en.W.values, grid, k),
                    b_up = aux_up[i].buoy[k],
                    b_en = en.B.values[k],
                    tke = prog_en.tke[k],
                    dMdz = ∇m,
                    M = m,
                    a_up = prog_up[i].area[k],
                    a_en = en.Area.values[k],
                    R_up = edmf.pressure_plume_spacing[i],
                    RH_up = aux_up[i].RH[k],
                    RH_en = en.RH.values[k],
                )

                er = entr_detr(param_set, εδ_model)
                edmf.entr_sc[i, k] = er.ε_dyn
                edmf.detr_sc[i, k] = er.δ_dyn
                # stochastic closure
                sde_model = edmf.sde_model
                stoch_ε = stochastic_closure(param_set, sde_model, Entrainment())
                stoch_δ = stochastic_closure(param_set, sde_model, Detrainment())
                edmf.entr_sc[i, k] *= stoch_ε
                edmf.detr_sc[i, k] *= stoch_δ

                edmf.frac_turb_entr[i, k] = er.ε_turb
                edmf.horiz_K_eddy[i, k] = er.K_ε
            else
                edmf.entr_sc[i, k] = 0.0
                edmf.detr_sc[i, k] = 0.0
                edmf.frac_turb_entr[i, k] = 0.0
                edmf.horiz_K_eddy[i, k] = 0.0
            end
        end
    end

    @inbounds for k in real_face_indices(grid)
        @inbounds for i in 1:(up.n_updrafts)

            # pressure
            a_bcs = (; bottom = SetValue(edmf.area_surface_bc[i]), top = SetValue(0))
            a_kfull = interpc2f(prog_up[i].area, grid, k; a_bcs...)
            if a_kfull > 0.0
                B = aux_up[i].buoy
                b_bcs = (; bottom = SetValue(B[kc_surf]), top = SetValue(B[kc_toa]))
                b_kfull = interpc2f(aux_up[i].buoy, grid, k; b_bcs...)
                w_cut = fcut(prog_up_f[i].w, grid, k)
                ∇w_up = f∇(w_cut, grid, k; bottom = SetValue(0), top = SetGradient(0))
                asp_ratio = 1.0
                nh_pressure_b, nh_pressure_adv, nh_pressure_drag = perturbation_pressure(
                    param_set,
                    up.updraft_top[i],
                    a_kfull,
                    b_kfull,
                    ρ0_f[k],
                    prog_up_f[i].w[k],
                    ∇w_up,
                    en.W.values[k],
                    asp_ratio,
                )
            else
                nh_pressure_b = 0.0
                nh_pressure_adv = 0.0
                nh_pressure_drag = 0.0
            end
            edmf.nh_pressure_b[i, k] = nh_pressure_b
            edmf.nh_pressure_adv[i, k] = nh_pressure_adv
            edmf.nh_pressure_drag[i, k] = nh_pressure_drag
            edmf.nh_pressure[i, k] = nh_pressure_b + nh_pressure_adv + nh_pressure_drag
        end
    end

    #####
    ##### compute_eddy_diffusivities_tke
    #####

    @inbounds for k in real_center_indices(grid)

        # compute shear
        U_cut = ccut(prog_gm.u, grid, k)
        V_cut = ccut(prog_gm.v, grid, k)
        wc_en = interpf2c(en.W.values, grid, k)
        wc_up = ntuple(up.n_updrafts) do i
            interpf2c(prog_up_f[i].w, grid, k)
        end
        w_dual = dual_faces(en.W.values, grid, k)

        ∇U = c∇(U_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇V = c∇(V_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        ∇w = ∇f2c(w_dual, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        Shear² = ∇U^2 + ∇V^2 + ∇w^2

        QT_cut = ccut(en.QT.values, grid, k)
        ∂qt∂z = c∇(QT_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        THL_cut = ccut(en.H.values, grid, k)
        ∂θl∂z = c∇(THL_cut, grid, k; bottom = SetGradient(0), top = SetGradient(0))
        # buoyancy_gradients
        bg_model = Tan2018(;
            qt_dry = en_thermo.qt_dry[k],
            th_dry = en_thermo.th_dry[k],
            t_cloudy = en_thermo.t_cloudy[k],
            qv_cloudy = en_thermo.qv_cloudy[k],
            qt_cloudy = en_thermo.qt_cloudy[k],
            th_cloudy = en_thermo.th_cloudy[k],
            ∂qt∂z = ∂qt∂z,
            ∂θl∂z = ∂θl∂z,
            p0 = p0_c[k],
            en_cld_frac = en.cloud_fraction.values[k],
            alpha0 = α0_c[k],
        )
        bg = buoyancy_gradients(param_set, bg_model)

        # Limiting stratification scale (Deardorff, 1976)
        p0_cut = ccut(p0_c, grid, k)
        T_cut = ccut(en.T.values, grid, k)
        QT_cut = ccut(en.QT.values, grid, k)
        QL_cut = ccut(en.QL.values, grid, k)
        ts_cut = TD.PhaseEquil_pTq.(param_set, p0_cut, T_cut, QT_cut)
        thv_cut = TD.virtual_pottemp.(ts_cut)

        ts = TD.PhaseEquil_pθq(param_set, p0_c[k], en.H.values[k], en.QT.values[k])
        θv = TD.virtual_pottemp(ts)
        ∂θv∂z = c∇(thv_cut, grid, k; bottom = SetGradient(0), top = Extrapolate())
        # compute ∇Ri and Pr
        ∇_Ri = gradient_Richardson_number(bg.∂b∂z_θl, bg.∂b∂z_qt, Shear², eps(0.0))
        edmf.prandtl_nvec[k] = turbulent_Prandtl_number(param_set, obukhov_length, ∇_Ri)

        ml_model = MinDisspLen(;
            z = grid.zc[k].z,
            obukhov_length = obukhov_length,
            tke_surf = prog_en.tke[kc_surf],
            ustar = surface.ustar,
            Pr = edmf.prandtl_nvec[k],
            p0 = p0_c[k],
            ∂b∂z_θl = bg.∂b∂z_θl,
            Shear² = Shear²,
            ∂b∂z_qt = bg.∂b∂z_qt,
            ∂θv∂z = ∂θv∂z,
            ∂qt∂z = ∂qt∂z,
            ∂θl∂z = ∂θl∂z,
            θv = θv,
            tke = prog_en.tke[k],
            a_en = (1 - aux_tc.bulk.area[k]),
            wc_en = wc_en,
            wc_up = Tuple(wc_up),
            a_up = ntuple(i -> prog_up[i].area[k], up.n_updrafts),
            ε_turb = ntuple(i -> edmf.frac_turb_entr[i, k], up.n_updrafts),
            δ_dyn = ntuple(i -> edmf.detr_sc[i, k], up.n_updrafts),
            en_cld_frac = en.cloud_fraction.values[k],
            θ_li_en = en.H.values[k],
            ql_en = en.QL.values[k],
            qt_en = en.QT.values[k],
            T_en = en.T.values[k],
            N_up = up.n_updrafts,
        )

        ml = mixing_length(param_set, ml_model)
        edmf.mls[k] = ml.min_len_ind
        edmf.mixing_length[k] = ml.mixing_length
        edmf.ml_ratio[k] = ml.ml_ratio

        KM[k] = c_m * edmf.mixing_length[k] * sqrt(max(prog_en.tke[k], 0.0))
        KH[k] = KM[k] / edmf.prandtl_nvec[k]

        aux_en_2m.tke.buoy[k] = -ml_model.a_en * ρ0_c[k] * KH[k] * (bg.∂b∂z_θl + bg.∂b∂z_qt)
    end

    compute_covariance_entr(edmf, grid, state, :tke, en.W, en.W, :w)
    compute_covariance_shear(edmf, grid, state, gm, :tke, en.W.values, en.W.values)
    compute_covariance_interdomain_src(edmf, grid, state, en.W, en.W, :tke, :w)

    #####
    ##### compute_tke_pressure
    #####
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.press[k] = 0.0
        @inbounds for i in 1:(up.n_updrafts)
            w_up_c = interpf2c(prog_up_f[i].w, grid, k)
            w_en_c = interpf2c(en.W.values, grid, k)
            press_c = interpf2c(edmf.nh_pressure, grid, k, i)
            aux_en_2m.tke.press[k] += (w_en_c - w_up_c) * press_c
        end
    end

    compute_covariance_entr(edmf, grid, state, :Hvar, en.H, en.H, :θ_liq_ice)
    compute_covariance_entr(edmf, grid, state, :QTvar, en.QT, en.QT, :q_tot)
    compute_covariance_entr(edmf, grid, state, :HQTcov, en.H, en.QT, :θ_liq_ice, :q_tot)
    compute_covariance_shear(edmf, grid, state, gm, :Hvar, en.H.values, en.H.values)
    compute_covariance_shear(edmf, grid, state, gm, :QTvar, en.QT.values, en.QT.values)
    compute_covariance_shear(edmf, grid, state, gm, :HQTcov, en.H.values, en.QT.values)
    compute_covariance_interdomain_src(edmf, grid, state, en.H, en.H, :Hvar, :θ_liq_ice)
    compute_covariance_interdomain_src(edmf, grid, state, en.QT, en.QT, :QTvar, :q_tot)
    compute_covariance_interdomain_src(edmf, grid, state, en.H, en.QT, :HQTcov, :θ_liq_ice, :q_tot)

    #####
    ##### compute_covariance_rain # need to update this one
    #####

    # TODO defined again in compute_covariance_shear and compute_covaraince
    ae = 1 .- aux_tc.bulk.area # area of environment
    @inbounds for k in real_center_indices(grid)
        aux_en_2m.tke.rain_src[k] = 0
        aux_en_2m.Hvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.Hvar_rain_dt[k]
        aux_en_2m.QTvar.rain_src[k] = ρ0_c[k] * ae[k] * 2 * en_thermo.QTvar_rain_dt[k]
        aux_en_2m.HQTcov.rain_src[k] = ρ0_c[k] * ae[k] * en_thermo.HQTcov_rain_dt[k]
    end

    reset_surface_covariance(edmf, grid, state, gm, Case)

    compute_diffusive_fluxes(edmf, grid, state, gm, Case, TS, param_set)
    update_cloud_frac(edmf, grid, state, gm)
end

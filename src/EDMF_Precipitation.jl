"""
Computes the qr advection (down) tendency
"""
function compute_rain_advection_tendencies(::PrecipPhysics, grid, state, gm, TS::TimeStepping)
    param_set = parameter_set(gm)
    Δz = grid.Δz
    Δt = TS.dt
    CFL_limit = 0.5

    ρ_0_c = center_ref_state(state).ρ0
    tendencies_pr = center_tendencies_precipitation(state)
    prog_pr = center_prog_precipitation(state)
    aux_tc = center_aux_turbconv(state)

    # helper to calculate the rain velocity
    # TODO: assuming gm.W = 0
    # TODO: verify translation
    term_vel_rain = aux_tc.term_vel_rain
    term_vel_snow = aux_tc.term_vel_snow

    @inbounds for k in real_center_indices(grid)
        term_vel_rain[k] = CM1.terminal_velocity(param_set, rain_type, ρ_0_c[k], prog_pr.qr[k])
        term_vel_snow[k] = CM1.terminal_velocity(param_set, snow_type, ρ_0_c[k], prog_pr.qs[k])
    end

    @inbounds for k in reverse(real_center_indices(grid))
        # check stability criterion
        CFL_out_rain = Δt / Δz * term_vel_rain[k]
        CFL_out_snow = Δt / Δz * term_vel_snow[k]
        if is_toa_center(grid, k)
            CFL_in_rain = 0.0
            CFL_in_snow = 0.0
        else
            CFL_in_rain = Δt / Δz * term_vel_rain[k + 1]
            CFL_in_snow = Δt / Δz * term_vel_snow[k + 1]
        end
        if max(CFL_in_rain, CFL_in_snow, CFL_out_rain, CFL_out_snow) > CFL_limit
            error("Time step is too large for rain fall velocity!")
        end

        ρ_0_cut = ccut_downwind(ρ_0_c, grid, k)

        QR_cut = ccut_downwind(prog_pr.qr, grid, k)
        QS_cut = ccut_downwind(prog_pr.qs, grid, k)

        w_cut_rain = ccut_downwind(term_vel_rain, grid, k)
        w_cut_snow = ccut_downwind(term_vel_snow, grid, k)

        ρQRw_cut = ρ_0_cut .* QR_cut .* w_cut_rain
        ρQSw_cut = ρ_0_cut .* QS_cut .* w_cut_snow

        ∇ρQRw = c∇_downwind(ρQRw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))
        ∇ρQSw = c∇_downwind(ρQSw_cut, grid, k; bottom = FreeBoundary(), top = SetValue(0))

        ρ_0_c_k = ρ_0_c[k]

        # TODO - some positivity limiters are needed
        aux_tc.qr_tendency_advection[k] = ∇ρQRw / ρ_0_c_k
        aux_tc.qs_tendency_advection[k] = ∇ρQSw / ρ_0_c_k

        tendencies_pr.qr[k] += aux_tc.qr_tendency_advection[k]
        tendencies_pr.qs[k] += aux_tc.qs_tendency_advection[k]
    end
    return
end

"""
Computes the tendencies to θ_liq_ice, qt and qr due to rain evaporation
"""
function compute_rain_evap_tendencies(::PrecipPhysics, grid, state, gm, TS::TimeStepping)
    param_set = parameter_set(gm)
    Δt = TS.dt
    p0_c = center_ref_state(state).p0
    ρ0_c = center_ref_state(state).ρ0
    aux_gm = center_aux_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    prog_gm = center_prog_grid_mean(state)
    prog_pr = center_prog_precipitation(state)
    tendencies_pr = center_tendencies_precipitation(state)

    @inbounds for k in real_center_indices(grid)
        q_tot_gm = prog_gm.q_tot[k]
        T_gm = aux_gm.T[k]
        # When we fuse loops, this should hopefully disappear
        ts = TD.PhaseEquil_pTq(param_set, p0_c[k], T_gm, q_tot_gm)
        q = TD.PhasePartition(ts)

        # TODO - move limiters elsewhere
        qt_tendency =
            min(prog_pr.qr[k] / Δt, -CM1.evaporation_sublimation(param_set, rain_type, q, prog_pr.qr[k], ρ0_c[k], T_gm))

        aux_tc.qt_tendency_rain_evap[k] = qt_tendency
        aux_tc.qr_tendency_rain_evap[k] = -qt_tendency
        tendencies_pr.qr[k] += aux_tc.qr_tendency_rain_evap[k]
        aux_tc.θ_liq_ice_tendency_rain_evap[k] = θ_liq_ice_helper(ts, qt_tendency)
    end
    return
end

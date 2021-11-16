#####
##### BaseCase methods
#####

function free_convection_windspeed(surf::SurfaceBase, grid, state, gm::GridMeanVariables, param_set, ::BaseCase)
    prog_gm = center_prog_grid_mean(state)
    aux_tc = center_aux_turbconv(state)
    zi = get_inversion(grid, state, param_set, surf.Ri_bulk_crit)
    wstar = get_wstar(surf.bflux, zi) # yair here zi in TRMM should be adjusted
    surf.windspeed = sqrt(surf.windspeed * surf.windspeed + (wstar) * (wstar))
    return
end

#####
##### Default SurfaceBase behavior
#####

free_convection_windspeed(surf::SurfaceBase, grid, state, gm::GridMeanVariables, param_set) =
    free_convection_windspeed(surf, grid, state, gm, param_set, BaseCase())

#####
##### SurfaceNone
#####

function update(surf::SurfaceBase{SurfaceNone}, grid, state, gm::GridMeanVariables, param_set)
    # JH: assigning small fluxed so that simulation won"t crash when computing mixing length
    # UPDATE  - YC  - the windspeed, zrough, bflux values that were here amounted to a
    # prescribed ustar 0.000540989171658113.
    # ustar plays a role in the surface covarianve and TKE and expalin the changes seen when
    # switching from SCAMPy to TC.jl in the second moment profile in the buble case.
    # For now I will prescribe these value in the Case.Surface instead of computing it in every timestep
    kc_surf = kc_surface(grid)
    return
end
free_convection_windspeed(surf::SurfaceBase{SurfaceNone}, grid, state, gm::GridMeanVariables, param_set) = nothing

function update(surf::SurfaceBase{SurfaceFixedFlux}, grid, state, gm::GridMeanVariables, param_set)
    FT = eltype(grid)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    ts_sfc = thermo_state_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)

    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts_in)
    zi = get_inversion(grid, state, param_set, surf.Ri_bulk_crit)
    gustiness = get_wstar(surf.bflux, zi) # yair here zi in TRMM should be adjusted

    u_sfc = SA.SVector{2,FT}(0, 0)
    u_in = SA.SVector{2,FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.ValuesSurface(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.ValuesInterior(z_in, u_in, ts_in)
    if surf.ustar_fixed
        sc = SF.FluxesAndFrictionVelocity{FT}(;
                state_in = vals_int,
                state_sfc = vals_sfc,
                shf = surf.shf,
                lhf = surf.lhf,
                ustar = surf.ustar,
                z0m = surf.zrough,
                z0b = surf.zrough,
                gustiness = gustiness,
            )
        result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
    else
        sc = SF.Fluxes{FT}(
                state_in = vals_int,
                state_sfc = vals_sfc,
                shf = surf.shf,
                lhf = surf.lhf,
                z0m = surf.zrough,
                z0b = surf.zrough,
                gustiness = gustiness,
            )
        result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
        surf.ustar = result.ustar
    end
    surf.obukhov_length = result.L_MO
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux =  surf.shf/TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf/TD.latent_heat_vapor(ts_in)
    return
end

function update(surf::SurfaceBase{SurfaceFixedCoeffs}, grid, state, gm::GridMeanVariables, param_set)
    FT = eltype(grid)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    aux_gm = center_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    T_gm_surf = aux_gm.T[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]

    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_sfc = thermo_state_pθq(param_set, p0_f_surf, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_f_surf, θ_liq_ice_gm_surf, q_tot_gm_surf)
    u_sfc = SA.SVector{2,FT}(0, 0)
    u_in = SA.SVector{2,FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.ValuesSurface(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.ValuesInterior(z_in, u_in, ts_in)
    sc = SF.Coefficients{FT}(
            state_in = vals_int,
            state_sfc = vals_sfc,
            Cd = surf.cm,
            Ch = surf.ch,
            z0m = surf.zrough,
            z0b = surf.zrough,
        )
    result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux = result.ρτxz
    surf.rho_vflux = result.ρτyz
    surf.rho_hflux =  surf.shf/TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf/TD.latent_heat_vapor(ts_in)
    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts_in)
    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhov}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    T_gm_surf = aux_gm.T[kc_surf]
    Pg = surf.ref_params.Pg

    ts_sfc = thermo_state_pTq(param_set, Pg, surf.Tsurface, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)

    u_sfc = SA.SVector{2,FT}(0, 0)
    u_in = SA.SVector{2,FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.ValuesSurface(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.ValuesInterior(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(
            state_in = vals_int,
            state_sfc = vals_sfc,
            z0m = surf.zrough,
            z0b = surf.zrough,
        )
    result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf * 0.0
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux  = result.ρτxz
    surf.rho_vflux  = result.ρτyz
    surf.rho_hflux =  surf.shf/TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf/TD.latent_heat_vapor(ts_in)
    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts_in)
    return
end

function update(surf::SurfaceBase{SurfaceMoninObukhovDry}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    T_gm_surf = aux_gm.T[kc_surf]
    Pg = surf.ref_params.Pg

    phase_part = TD.PhasePartition(surf.qsurface, 0.0, 0.0)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, Pg, phase_part)
    ts_sfc = thermo_state_pθq(param_set, Pg, h_star, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)
    u_sfc = SA.SVector{2,FT}(0, 0)
    u_in = SA.SVector{2,FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.ValuesSurface(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.ValuesInterior(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(
            state_in = vals_int,
            state_sfc = vals_sfc,
            z0m = surf.zrough,
            z0b = surf.zrough,
        )
    result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf * 0.0
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux  = result.ρτxz
    surf.rho_vflux  = result.ρτyz
    surf.rho_hflux =  surf.shf/TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf/TD.latent_heat_vapor(ts_in)
    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts_in)
    return
end

function update(surf::SurfaceBase{SurfaceSullivanPatton}, grid, state, gm::GridMeanVariables, param_set)
    g = CPP.grav(param_set)
    kc_surf = kc_surface(grid)
    kf_surf = kf_surface(grid)
    FT = eltype(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    p0_f_surf = face_ref_state(state).p0[kf_surf]
    p0_c_surf = center_ref_state(state).p0[kc_surf]
    ρ0_f_surf = face_ref_state(state).ρ0[kf_surf]
    α0_f_surf = face_ref_state(state).α0[kf_surf]
    prog_gm = center_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    u_gm_surf = prog_gm.u[kc_surf]
    v_gm_surf = prog_gm.v[kc_surf]
    q_tot_gm_surf = prog_gm.q_tot[kc_surf]
    θ_liq_ice_gm_surf = prog_gm.θ_liq_ice[kc_surf]
    T_gm_surf = aux_gm.T[kc_surf]
    Pg = surf.ref_params.Pg

    phase_part = TD.PhasePartition(q_tot_gm_surf, 0.0, 0.0)
    pvg = TD.saturation_vapor_pressure(param_set, TD.PhaseEquil, surf.Tsurface)
    surf.qsurface = TD.q_vap_saturation_from_density(param_set, surf.Tsurface, ρ0_f_surf, pvg)
    h_star = TD.liquid_ice_pottemp_given_pressure(param_set, surf.Tsurface, Pg, phase_part)

    ts_sfc = thermo_state_pθq(param_set, Pg, h_star, surf.qsurface)
    ts_in = thermo_state_pθq(param_set, p0_c_surf, θ_liq_ice_gm_surf, surf.qsurface)

    u_sfc = SA.SVector{2,FT}(0, 0)
    u_in = SA.SVector{2,FT}(u_gm_surf, v_gm_surf)
    vals_sfc = SF.ValuesSurface(z_sfc, u_sfc, ts_sfc)
    vals_int = SF.ValuesInterior(z_in, u_in, ts_in)
    sc = SF.ValuesOnly{FT}(
            state_in = vals_int,
            state_sfc = vals_sfc,
            z0m = surf.zrough,
            z0b = surf.zrough,
        )
    result = SF.surface_conditions(param_set, sc, UF.Businger, SF.FDScheme())
    surf.cm = result.Cd
    surf.ch = result.Ch
    surf.obukhov_length = result.L_MO
    surf.lhf = result.lhf * 0.0
    surf.shf = result.shf
    surf.ustar = result.ustar
    surf.rho_uflux  = result.ρτxz
    surf.rho_vflux  = result.ρτyz
    surf.rho_hflux =  surf.shf/TD.cp_m(ts_in)
    surf.rho_qtflux = surf.lhf/TD.latent_heat_vapor(ts_in)
    surf.bflux = buoyancy_flux(param_set, surf.shf, surf.lhf, surf.Tsurface, surf.qsurface, α0_f_surf, ts_in)
    return
end

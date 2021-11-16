
function buoyancy_flux(param_set, shf::FT, lhf, T_b, qt_b, α0_b, ts) where {FT}
    g = FT(CPP.grav(param_set))
    molmass_ratio = FT(CPP.molmass_ratio(param_set))
    lv = TD.latent_heat_vapor(param_set, T_b)
    cp_m = TD.cp_m(ts)
    return (g * α0_b / cp_m / T_b * (shf + (molmass_ratio - 1) * cp_m * T_b * lhf / lv))
end
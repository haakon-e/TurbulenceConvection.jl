function remove_precipitation_rain(
    param_set::APS,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    τ_precip_rain::FT = 2500.0
    qc_0_rain::FT = 0.5e-3 #1.5 2.8 #3.0 #3.2 #3.5

    return -max(0, (q.liq - qc_0_rain)) / τ_precip_rain
end
function remove_precipitation_snow(
    param_set::APS,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    τ_precip_snow::FT = 100.0
    qc_0_snow::FT = 1e-6

    return -max(0, (q.ice - qc_0_snow)) / τ_precip_snow
end

"""
Computes the tendencies to qt and θ_liq_ice due to precipitation formation
(autoconversion + accretion)
"""
function precipitation_formation(param_set::APS, precipitation_model, qr::FT, qs::FT, area, ρ0, dt, ts) where {FT}

    # TODO - when using adaptive timestepping we are limiting the source terms
    #        with the previous timestep dt
    qt_tendency = 0.0
    qr_tendency = 0.0
    qs_tendency = 0.0
    θ_liq_ice_tendency = 0.0

    if area > 0.0

        q = TD.PhasePartition(ts)

        Π_m = TD.exner(ts)
        c_pm = TD.cp_m(ts)
        L_v0 = CPP.LH_v0(param_set)
        L_s0 = CPP.LH_s0(param_set)

        if precipitation_model == "cutoff"
            qsat = TD.q_vap_saturation(ts)
            λ = TD.liquid_fraction(ts)

            #S_qt = -min((q.liq + q.ice) / dt, -CM0.remove_precipitation(param_set, q, qsat))
            S_qt = -min((q.liq + q.ice) / dt, -CM0.remove_precipitation(param_set, q))

            qr_tendency -= S_qt * λ
            qs_tendency -= S_qt * (1.0 - λ)
            qt_tendency += S_qt
            θ_liq_ice_tendency -= S_qt / Π_m / c_pm * (L_v0 * λ + L_s0 * (1 - λ))
        end

        if precipitation_model == "clima_1m"
            T = TD.air_temperature(ts)
            T_fr = CPP.T_freeze(param_set)
            c_vl = CPP.cv_l(param_set)
            c_vm = TD.cv_m(ts)
            Rm = TD.gas_constant_air(ts)
            Lf = TD.latent_heat_fusion(ts)

            qsat = TD.q_vap_saturation(ts)
            λ = TD.liquid_fraction(ts)

            # autoconversion of cloud water and ice handled using the 0-moment scheme rates
            # saturation adjustment prevents using the 1-moment snow autoconversion rate
            # because it assumes supersaturation is present
            #S_qt = -min((q.liq + q.ice) / dt, -CM0.remove_precipitation(param_set, q, qsat))
            S_qt_rain = -min(q.liq / dt, -remove_precipitation_rain(param_set, q))
            S_qt_snow = -min(q.ice / dt, -remove_precipitation_snow(param_set, q))

            qr_tendency -= S_qt_rain
            qs_tendency -= S_qt_snow
            qt_tendency += S_qt_rain + S_qt_snow
            θ_liq_ice_tendency -= 1.0 / Π_m / c_pm * (L_v0 * S_qt_rain + L_s0 * S_qt_snow)

            # accretion cloud water + rain
            S_qr = min(q.liq / dt, CM1.accretion(param_set, liq_type, rain_type, q.liq, qr, ρ0))
            qr_tendency += S_qr
            qt_tendency -= S_qr
            θ_liq_ice_tendency += S_qr / Π_m / c_pm * L_v0

            # accretion cloud ice + snow
            S_qs = min(q.ice / dt, CM1.accretion(param_set, ice_type, snow_type, q.ice, qs, ρ0))
            qs_tendency += S_qs
            qt_tendency -= S_qs
            θ_liq_ice_tendency += S_qs / Π_m / c_pm * L_s0

            # sink of cloud water via accretion cloud water + snow
            S_qt = -min(q.liq / dt, CM1.accretion(param_set, liq_type, snow_type, q.liq, qs, ρ0))
            if T < T_fr # cloud droplets freeze to become snow)
                qs_tendency -= S_qt
                qt_tendency += S_qt
                θ_liq_ice_tendency -= S_qt / Π_m / c_pm * Lf * (1.0 + Rm / c_vm)
            else # snow melts, both cloud water and snow become rain
                α::FT = c_vl / Lf * (T - T_fr)
                qt_tendency += S_qt
                qs_tendency += S_qt * α
                qr_tendency -= S_qt * (1 + α)
                θ_liq_ice_tendency += S_qt / Π_m / c_pm * (Lf * (1.0 + Rm / c_vm) * α - L_v0)
            end

            # sink of cloud ice via accretion cloud ice - rain
            S_qt = -min(q.ice / dt, CM1.accretion(param_set, ice_type, rain_type, q.ice, qr, ρ0))
            # sink of rain via accretion cloud ice - rain
            S_qr = -min(qr / dt, CM1.accretion_rain_sink(param_set, q.ice, qr, ρ0))
            qt_tendency += S_qt
            qr_tendency += S_qr
            qs_tendency += -(S_qt + S_qr)
            θ_liq_ice_tendency -= 1.0 / Π_m / c_pm * (S_qr * Lf * (1.0 + Rm / c_vm) + S_qt * L_s0)

            # accretion rain - snow
            if T < T_fr
                S_qs = min(qr / dt, CM1.accretion_snow_rain(param_set, snow_type, rain_type, qs, qr, ρ0))
            else
                S_qs = -min(qs / dt, CM1.accretion_snow_rain(param_set, rain_type, snow_type, qr, qs, ρ0))
            end
            qs_tendency += S_qs
            qr_tendency -= S_qs
            θ_liq_ice_tendency += S_qs * Lf / Π_m / c_vm
        end
    end
    return PrecipFormation(θ_liq_ice_tendency, qt_tendency, qr_tendency, qs_tendency)
end

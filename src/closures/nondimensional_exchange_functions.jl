"""
    nondimensional_exchange_functions(
        param_set,
        c_δ,
        c_μ,
        c_μ0,
        β,
        χ_up,
        Δw,
        Δb,
        εδ_model,
    )

Returns the nondimensional entrainment and detrainment
functions following Cohen et al. (JAMES, 2020), given:
 - `param_set` a parameter set
 - `c_δ`,  detrainment factor
 - `c_μ`,  logistic function factor
 - `c_μ0`, logistic function timescale (sec)
 - `β`,    sorting power
 - `χ_up`, updraft mixing fraction
 - `Δw`, updraft - environment vertical velocity differnce
 - `Δb`, updraft - environment buoynacy differnce
 - `εδ_model`, entrainment detrainment model type
"""
function nondimensional_exchange_functions(param_set, c_δ, c_μ, c_μ0, β, χ_upd, Δw, Δb, εδ_model::MoistureDeficitEntr)

    # should be: c_δ = sign(condensate(ts_en) + condensate(ts_up[i])) * entr.c_δ
    μ_ij = (χ_upd - εδ_model.a_up / (εδ_model.a_up + εδ_model.a_en)) * Δb / Δw
    D_ε = 1.0 / (1.0 + exp(-c_μ / c_μ0 * μ_ij))
    D_δ = 1.0 / (1.0 + exp(c_μ / c_μ0 * μ_ij))

    M_δ = (max((εδ_model.RH_up / 100.0)^β - (εδ_model.RH_en / 100.0)^β, 0.0))^(1.0 / β)
    M_ε = (max((εδ_model.RH_en / 100.0)^β - (εδ_model.RH_up / 100.0)^β, 0.0))^(1.0 / β)

    return D_ε, D_δ, M_δ, M_ε
end;

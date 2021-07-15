#! format: off
var_map_les(s::String) = var_map_les(Val(Symbol(s)))
var_map_les(s::Symbol) = var_map_les(Val(s))
var_map_les(::Val{T}) where {T} = nothing

var_map_scampy(s::String) = s

# var_map_les(::Val{:tc_var}) = "les_var"
var_map_les(::Val{:rho}) = "rho"
var_map_les(::Val{:u_mean}) = "u_mean"
var_map_les(::Val{:v_mean}) = "v_mean"
var_map_les(::Val{:qt_mean}) = "qt_mean"
var_map_les(::Val{:ql_mean}) = "ql_mean"
# var_map_les(::Val{:updraft_fraction}) = "updraft_area"
var_map_les(::Val{:updraft_area}) = "updraft_fraction"
var_map_les(::Val{:updraft_w}) = "updraft_w"
var_map_les(::Val{:updraft_qt}) = "updraft_qt"
var_map_les(::Val{:updraft_ql}) = "updraft_ql"
var_map_les(::Val{:updraft_qr}) = "updraft_ql"
var_map_les(::Val{:updraft_thetal}) = "updraft_thetali"
var_map_les(::Val{:tke_mean}) = "tke_mean"
#! format: on
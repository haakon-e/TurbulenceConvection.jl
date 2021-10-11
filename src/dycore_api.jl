#####
##### Dycore API
#####

#=
This file provides a list of methods that TurbulenceConvection.jl
expects that the host dycore supports. This is experimental, as
we're not sure how the data structures / flow control will shake out.
=#

""" The cell center fields of any state vector """
center_prog(state) = state.cent

""" The cell face fields of any state vector """
face_prog(state) = state.face

""" The cell center reference state fields """
center_ref_state(state) = state.aux.cent.ref_state

""" The cell face reference state fields """
face_ref_state(state) = state.aux.face.ref_state

# Launch with `julia --project --track-allocation=user`
if !("." in LOAD_PATH)
    push!(LOAD_PATH, ".")
end
import Profile

include("common.jl")
case_name = ENV["ALLOCATION_CASE_NAME"]
@info "Recording allocations for $case_name"
sim = init_sim(case_name)
tendencies = copy(sim.state.prog)
main(NameList.default_namelist(case_name)) # compile first
Profile.clear_malloc_data()
main(NameList.default_namelist(case_name))

# Quit julia (which generates .mem files), then call
#=
import Coverage
allocs = Coverage.analyze_malloc("src")
=#

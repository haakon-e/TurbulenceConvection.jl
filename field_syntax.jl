import ClimaCore
const CC = ClimaCore
const CCO = ClimaCore.Operators

struct Grid{FT, CS, FS, SC, SF}
    zmin::FT
    zmax::FT
    Δz::FT
    Δzi::FT
    nz::Int
    cs::CS
    fs::FS
    zc::SC
    zf::SF
    function Grid(Δz::FT, nz::Int) where {FT <: AbstractFloat}
        z₀, z₁ = FT(0), FT(nz * Δz)
        domain = CC.Domains.IntervalDomain(CC.Geometry.ZPoint{FT}(z₀), CC.Geometry.ZPoint{FT}(z₁), boundary_tags = (:bottom, :top))
        mesh = CC.Meshes.IntervalMesh(domain, nelems = nz)
        cs = CC.Spaces.CenterFiniteDifferenceSpace(mesh); fs = CC.Spaces.FaceFiniteDifferenceSpace(cs)
        zc = CC.Fields.coordinate_field(cs); zf = CC.Fields.coordinate_field(fs)
        Δzi = 1 / Δz
        zmin = minimum(parent(zf)); zmax = maximum(parent(zf))
        CS = typeof(cs); FS = typeof(fs); SC = typeof(zc); SF = typeof(zf)
        return new{FT, CS, FS, SC, SF}(zmin, zmax, Δz, Δzi, nz, cs, fs, zc, zf)
    end
end
function FieldFromNamedTuple(space, nt::NamedTuple)
    cmv(z) = nt; return cmv.(CC.Fields.coordinate_field(space))
end

kc_surface(grid::Grid) = 1
kf_surface(grid::Grid) = 1
kc_top_of_atmos(grid::Grid) = grid.nz
kf_top_of_atmos(grid::Grid) = grid.nz + 1

real_center_indices(grid::Grid) = kc_surface(grid):kc_top_of_atmos(grid)
real_face_indices(grid::Grid) = kf_surface(grid):kf_top_of_atmos(grid)

prognostic_vars(FT, n_up) = (; prognostic_vars_gm(FT)..., prognostic_vars_edmf(FT, n_up)...)
prognostic_vars_gm(FT) = (; U = FT(0), V = FT(0), H = FT(0), QT = FT(0))
prognostic_vars_up(FT) = (; Area = FT(0), H = FT(0), QT = FT(0))
prognostic_vars_en(FT) = (; TKE = FT(0), Hvar = FT(0), QTvar = FT(0), HQTcov = FT(0))
prognostic_vars_edmf(FT, n_up) =
    (; turbconv = (; en = prognostic_vars_en(FT), up = ntuple(i -> prognostic_vars_up(FT), n_up)))
center_space(grid::Grid) = grid.cs
face_space(grid::Grid) = grid.fs

FT = Float64
n_updrafts = 2
grid = Grid(FT(0.1), 10)
cent_state_fields = FieldFromNamedTuple(center_space(grid), prognostic_vars(FT, n_updrafts))
face_state_fields = FieldFromNamedTuple(face_space(grid), prognostic_vars(FT, n_updrafts))
state = CC.Fields.FieldVector(cent = cent_state_fields, face = face_state_fields)

struct Cent{I <: Integer}
    i::I
end

Base.@propagate_inbounds Base.getindex(field::CC.Fields.FiniteDifferenceField, i::Integer) =
    Base.getproperty(field, i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.CenterFiniteDifferenceField, i::Cent) =
    Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.CenterFiniteDifferenceField, v, i::Cent) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.FaceFiniteDifferenceField, i::CCO.PlusHalf) =
    Base.getindex(CC.Fields.field_values(field), i.i)
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.FaceFiniteDifferenceField, v, i::CCO.PlusHalf) =
    Base.setindex!(CC.Fields.field_values(field), v, i.i)

Base.@propagate_inbounds Base.getindex(field::CC.Fields.FaceFiniteDifferenceField, ::Cent) =
    error("Attempting to getindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.getindex(field::CC.Fields.CenterFiniteDifferenceField, ::CCO.PlusHalf) =
    error("Attempting to getindex with a face index (PlusHalf) into a Center field")

Base.@propagate_inbounds Base.setindex!(field::CC.Fields.FaceFiniteDifferenceField, v, ::Cent) =
    error("Attempting to setindex with a center index (Cent) into a Face field")
Base.@propagate_inbounds Base.setindex!(field::CC.Fields.CenterFiniteDifferenceField, v, ::CCO.PlusHalf) =
    error("Attempting to setindex with a face index (PlusHalf) into a Center field")

using Test
for k in real_center_indices(grid)
    for i in 1:n_updrafts
        x = state.cent.turbconv.up[i].Area[Cent(k)]
        state.cent.turbconv.up[i].Area[Cent(k)] = 2

        @test_throws ErrorException x = state.face.turbconv.up[i].Area[Cent(k)]
        @test_throws ErrorException state.face.turbconv.up[i].Area[Cent(k)] = 2

        @test_throws ErrorException x = state.cent.turbconv.up[i].Area[CCO.PlusHalf(k)]
        @test_throws ErrorException state.cent.turbconv.up[i].Area[CCO.PlusHalf(k)] = 2

        x = state.face.turbconv.up[i].Area[CCO.PlusHalf(k)]
        state.face.turbconv.up[i].Area[CCO.PlusHalf(k)] = 2
        state.face.turbconv.up[i].QT[CCO.PlusHalf(k)] = 3
    end
end

# state.cent.turbconv.up[i].Area[k] = x
# state.cent.turbconv.up[:].Area[k] = ntuple(...)

state.cent
state.face



# TODO: remove `vars` hack that avoids https://github.com/Alexander-Barth/NCDatasets.jl/issues/135

mutable struct NetCDFIO_Stats
    root_grp::NC.NCDataset{Nothing}
    profiles_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    ts_grp::NC.NCDataset{NC.NCDataset{Nothing}}
    grid::Grid
    last_output_time::Float64
    uuid::String
    frequency::Float64
    stats_path::String
    path_plus_file::String
    vars::Dict{String, Any} # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    function NetCDFIO_Stats(namelist, grid::Grid)

        # Initialize properties with valid type:
        tmp = tempname()
        root_grp = NC.Dataset(tmp, "c")
        NC.defGroup(root_grp, "profiles")
        NC.defGroup(root_grp, "timeseries")
        profiles_grp = root_grp.group["profiles"]
        ts_grp = root_grp.group["timeseries"]
        close(root_grp)

        last_output_time = 0.0
        uuid = string(namelist["meta"]["uuid"])

        frequency = namelist["stats_io"]["frequency"]

        # Setup the statistics output path
        simname = namelist["meta"]["simname"]
        casename = namelist["meta"]["casename"]
        outpath = joinpath(namelist["output"]["output_root"], "Output.$simname.$uuid")
        mkpath(outpath)

        stats_path = joinpath(outpath, namelist["stats_io"]["stats_dir"])
        mkpath(stats_path)

        path_plus_file = joinpath(stats_path, "Stats.$simname.nc")

        # TODO: uncomment restart
        # if isfile(path_plus_file)
        #   @inbounds for i in 1:100
        #         res_name = "Restart_$i"
        #         if isfile(path_plus_file)
        #             path_plus_file = stats_path * "Stats.$simname.$res_name.nc"
        #         else
        #             break
        #         end
        #     end
        # end

        # Write namelist file to output directory
        open(joinpath(outpath, "namelist_$casename.in"), "w") do io
            JSON.print(io, namelist, 4)
        end

        # Remove the NC file if it exists, in case it accidentally wasn't closed
        isfile(path_plus_file) && rm(path_plus_file; force = true)

        NC.Dataset(path_plus_file, "c") do root_grp

            zf = vec(grid.zf)
            zc = vec(grid.zc)

            # Set profile dimensions
            profile_grp = NC.defGroup(root_grp, "profiles")
            NC.defDim(profile_grp, "zf", grid.nz + 1)
            NC.defDim(profile_grp, "zc", grid.nz)
            NC.defDim(profile_grp, "t", Inf)
            NC.defVar(profile_grp, "zf", zf, ("zf",))
            NC.defVar(profile_grp, "zc", zc, ("zc",))
            NC.defVar(profile_grp, "t", Float64, ("t",))

            reference_grp = NC.defGroup(root_grp, "reference")
            NC.defDim(reference_grp, "zf", grid.nz + 1)
            NC.defDim(reference_grp, "zc", grid.nz)
            NC.defVar(reference_grp, "zf", zf, ("zf",))
            NC.defVar(reference_grp, "zc", zc, ("zc",))

            ts_grp = NC.defGroup(root_grp, "timeseries")
            NC.defDim(ts_grp, "t", Inf)
            NC.defVar(ts_grp, "t", Float64, ("t",))
        end
        vars = Dict{String, Any}()
        return new(
            root_grp,
            profiles_grp,
            ts_grp,
            grid,
            last_output_time,
            uuid,
            frequency,
            stats_path,
            path_plus_file,
            vars,
        )
    end
end


function open_files(self)
    self.root_grp = NC.Dataset(self.path_plus_file, "a")
    self.profiles_grp = self.root_grp.group["profiles"]
    self.ts_grp = self.root_grp.group["timeseries"]
    vars = self.vars

    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    vars["profiles"] = Dict{String, Any}()
    for k in keys(self.profiles_grp)
        vars["profiles"][k] = self.profiles_grp[k]
    end
    vars["timeseries"] = Dict{String, Any}()
    for k in keys(self.ts_grp)
        vars["timeseries"][k] = self.ts_grp[k]
    end
end

function close_files(self::NetCDFIO_Stats)
    close(self.root_grp)
end

#####
##### Generic field
#####

function add_field(self::NetCDFIO_Stats, var_name::String; dims, group)
    NC.Dataset(self.path_plus_file, "a") do root_grp
        profile_grp = root_grp.group[group]
        new_var = NC.defVar(profile_grp, var_name, Float64, dims)
    end
end

#####
##### Time-series data
#####

function add_ts(self::NetCDFIO_Stats, var_name::String)
    NC.Dataset(self.path_plus_file, "a") do root_grp
        ts_grp = root_grp.group["timeseries"]
        new_var = NC.defVar(ts_grp, var_name, Float64, ("t",))
    end
end

#####
##### Performance critical IO
#####

# Field wrapper
write_field(self::NetCDFIO_Stats, var_name::String, data; group) = write_field(self, var_name, vec(data); group = group)

function write_field(self::NetCDFIO_Stats, var_name::String, data::T; group) where {T <: AbstractArray{Float64, 1}}
    if group == "profiles"
        # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
        @inbounds self.vars[group][var_name][:, end] = data
        # Ideally, we remove self.vars and use:
        # var = self.profiles_grp[var_name]
        # Not sure why `end` instead of `end+1`, but `end+1` produces garbage output
        # @inbounds var[end, :] = data :: T
    elseif group == "reference"
        NC.Dataset(self.path_plus_file, "a") do root_grp
            reference_grp = root_grp.group[group]
            var = reference_grp[var_name]
            var .= data::T
        end
    else
        error("Bad group given")
    end
end

function write_ts(self::NetCDFIO_Stats, var_name::String, data::Float64)
    # Hack to avoid https://github.com/Alexander-Barth/NCDatasets.jl/issues/135
    @inbounds self.vars["timeseries"][var_name][end] = data::Float64
    # Ideally, we remove self.vars and use:
    # var = self.ts_grp[var_name]
    # @inbounds var[end+1] = data :: Float64
end

function write_simulation_time(self::NetCDFIO_Stats, t::Float64)
    # # Write to profiles group
    profile_t = self.profiles_grp["t"]
    @inbounds profile_t[end + 1] = t::Float64

    # # Write to timeseries group
    ts_t = self.ts_grp["t"]
    @inbounds ts_t[end + 1] = t::Float64
end

module ClimaSeaIceNCDatasetsExt

using NCDatasets
using Oceananigans
using ClimaSeaIce

#####
##### Variable attributes
#####

default_horizontal_velocity_attributes(::RectilinearGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the +x-direction.", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the +y-direction.", "units" => "m/s"))

default_horizontal_velocity_attributes(::LatitudeLongitudeGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the zonal direction (+ = east).", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the meridional direction (+ = north).", "units" => "m/s"))

default_horizontal_velocity_attributes(::OrthogonalSphericalShellGrid) = Dict(
    "u" => Dict("long_name" => "Velocity in the i-direction (+ = increasing i).", "units" => "m/s"),
    "v" => Dict("long_name" => "Velocity in the j-direction (+ = increasing j).", "units" => "m/s")
)
    
default_horizontal_velocity_attributes(ibg::ImmersedBoundaryGrid) = default_horizontal_velocity_attributes(ibg.underlying_grid)

default_sea_ice_attributes() = Dict(
    "h"  => Dict("long_name" => "Sea ice thickness.", "units" => "m"),
    "ℵ"  => Dict("long_name" => "Sea ice concentration.", "units" => "-"),
    "hs" => Dict("long_name" => "Snow thickness.", "units" => "m")
)

function __init__()
    OCNE = Base.get_extension(Oceananigans, :OceananigansNCDatasetsExt)
    if !isnothing(OCNE)
        @eval begin
            function $OCNE.default_output_attributes(model::SeaIceModel)
                velocity_attrs = default_horizontal_velocity_attributes(model.grid)
                tracer_attrs = default_sea_ice_attributes()
                return merge(velocity_attrs, tracer_attrs)
            end
        end
    end
end

end

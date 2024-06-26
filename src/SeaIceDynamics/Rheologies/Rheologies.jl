module Rheologies

using KernelAbstractions: @kernel, @index

using Adapt

using ClimaSeaIce.SeaIceDynamics:
    AbstractRheology,
    update_stepping_coefficients!,
    get_stepping_coefficients

include("nothing_rheology.jl")
include("explicit_visco_plastic_rheology.jl")

end

module ClimaSeaIce

include("EulerianThermodynamicSeaIceModels.jl")
include("LagrangianThermodynamicSeaIceModels.jl")

using .EulerianThermodynamicSeaIceModels: EulerianThermodynamicSeaIceModel

end # module

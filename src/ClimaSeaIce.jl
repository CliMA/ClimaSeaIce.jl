module ClimaSeaIce

include("EulerianThermodynamicSeaIceModels.jl")
include("PrognosticThicknessThermodynamicsModels.jl")

using .EulerianThermodynamicSeaIceModels: EulerianThermodynamicSeaIceModel

end # module

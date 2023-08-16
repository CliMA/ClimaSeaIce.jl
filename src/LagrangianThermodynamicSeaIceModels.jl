module LagrangianThermodynamicSeaIceModels

abstract type AbstractLagrangianThermodynamicSeaIceModel end

mutable struct ZeroLayerThermodynamicSeaIceModel{G, H, A, L, KS, KI} <: AbstractLagrangianThermodynamicSeaIceModel 
    horizontal_grid :: G
    thicknesses :: H
    previous_thicknesses :: H
    albedo :: A
    latent_heat_of_fusion :: L
    snow_diffusivity :: KS
    ice_diffusivity :: KI
end

end # module

module PrognosticThicknessThermodynamicsModels

abstract type AbstractPrognosticThicknessThermodynamicsModel end

mutable struct ZeroLayerThermodynamicSeaIceModel{G, H, A, L, KS, KI} <: AbstractPrognosticThicknessThermodynamicsModel
    horizontal_grid :: G
    thicknesses :: H
    previous_thicknesses :: H
    albedo :: A
    latent_heat_of_fusion :: L
    snow_diffusivity :: KS
    ice_diffusivity :: KI
end

end # module

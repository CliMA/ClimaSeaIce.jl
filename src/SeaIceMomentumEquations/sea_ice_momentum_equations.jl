
struct SeaIceMomentumEquation{C, R, T, B}
    coriolis :: C
    rheology :: R
    top_stress :: T
    bottom_stress :: B
end

function SeaIceMomentumEquation(; coriolis = nothing,
                                  rheology = nothing,
                                  top_stress = nothing,
                                  bottom_stress = nothing)

    return SeaIceMomentumEquation(coriolis, rheology, top_stress, bottom_stress)
end
    
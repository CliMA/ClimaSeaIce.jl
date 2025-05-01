
struct FluxBalanceFreeDrift{T, B, C}
    top_momentum_stess :: T
    bottom_momentum_stress :: B
end

# Just passing ocean velocities without mitigation
@inline function free_drift_u(i, j, k, grid, f::FluxBalanceFreeDrift, clock, fields)
    τib = implicit_τx_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τx_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τx(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τx(i, j, k, grid, f.top_momentum_stress, clock, fields)

    return (τib - τit) / (τeb - τet)
end

# Just passing ocean velocities without mitigation
@inline function free_drift_v(i, j, k, grid, f::FluxBalanceFreeDrift, clock, fields) 
    τib = implicit_τy_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τy_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τy(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τy(i, j, k, grid, f.top_momentum_stress, clock, fields)

    return (τib - τit) / (τeb - τet)
end
    

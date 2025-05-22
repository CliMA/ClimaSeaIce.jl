
"""
    StressBalanceFreeDrift{T, B, C}

A free drift parameterization that computes the free drift velocities as a balance between
top and bottom stresses ``τa ≈ τo`` where we split the stresses into a linear implicit part
and an explicit part  ``τaˣ = τaˣₑ + u * τaˣᵢ`` and ``τoˣ = τoˣₑ + u * τoˣᵢ`` (and similarly for the y-component) such that 
```
uᶠ = (τoˣₑ - τaˣₑ) / (τoˣᵢ - τaˣᵢ)
vᶠ = (τoʸₑ - τaʸₑ) / (τoʸᵢ - τaʸᵢ)
```
"""
struct StressBalanceFreeDrift{T, B, C}
    top_momentum_stess :: T
    bottom_momentum_stress :: B
end

# Just passing ocean velocities without mitigation
@inline function free_drift_u(i, j, k, grid, f::StressBalanceFreeDrift, clock, fields)
    τib = implicit_τx_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τx_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τx(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τx(i, j, k, grid, f.top_momentum_stress, clock, fields)

    # Explicit stress
    τe = τeb - τet

    # Implicit stress component 
    τi = τib - τit

    return ifelse(τe == 0, zero(grid), τi / τe)
end

# Just passing ocean velocities without mitigation
@inline function free_drift_v(i, j, k, grid, f::StressBalanceFreeDrift, clock, fields) 
    τib = implicit_τy_coefficient(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τit = implicit_τy_coefficient(i, j, k, grid, f.top_momentum_stress, clock, fields)

    τeb = explicit_τy(i, j, k, grid, f.bottom_momentum_stress, clock, fields)
    τet = explicit_τy(i, j, k, grid, f.top_momentum_stress, clock, fields)

    # Explicit stress
    τe = τeb - τet

    # Implicit stress component 
    τi = τib - τit

    return ifelse(τe == 0, zero(grid), τi / τe)
end
    

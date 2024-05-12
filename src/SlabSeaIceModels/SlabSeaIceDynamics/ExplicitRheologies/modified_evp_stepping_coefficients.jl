
struct ModifiedEVPSteppingCoefficients{A}
    c :: A
end

# Following the mEVP formulation of Kimmritz et al (2016)
function ModifiedEVPSteppingCoefficients(grid::AbstractGrid) 
    c = CenterField(grid)
    set!(c, 300)
    return ModifiedEVPSteppingCoefficients(c)
end

@inline get_stepping_coefficients(i, j, k, grid, rheology, coefficients) = rheology.substeps / 2
@inline get_stepping_coefficients(i, j, k, grid, rheology, coefficients::ModifiedEVPSteppingCoefficients) = @inbounds coefficients.c[i, j, k] 

@inline update_stepping_coefficients!(i, j, k, grid, coefficients, args...) = nothing

# Following an mEVP formulation: α = β = sqrt(4γ) (Kimmritz et al 2016)
# where γ = (ζ / 4) * π² * (Δt / mᵢ) / Az
@inline function update_stepping_coefficients!(i, j, k, grid, coefficients::ModifiedEVPSteppingCoefficients, ζ, mᵢ, Δt)
    A = Azᶜᶜᶜ(i, j, k, grid)
    γ = ζ / 4 * π^2 * Δt / mᵢ / A
    @inbounds coefficients.c[i, j, k] = sqrt(4γ)

    return nothing
end
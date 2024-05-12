struct ModifiedEVPSteppingCoefficients{A}
    c :: A
end

"""
    ModifiedEVPSteppingCoefficients(grid::AbstractGrid)

Constructs a `ModifiedEVPSteppingCoefficients` object based on the given grid that
represents the substepping coefficients for the modified EVP formulation of Kimmritz et al (2016).

In particular: α (stress stepping) == β (momentum stepping) = sqrt(4γ) 
where γ = (ζ / 4) * π² * (Δt / mᵢ) / Az is a stability parameter 

This formulation allows fast convergence in regions where sqrt(4γ) is small. Regions where
the coefficients are large correspond to regions where the ice is more solid and moves as a block
and the convergence is slower.
"""
function ModifiedEVPSteppingCoefficients(grid::AbstractGrid) 
    c = CenterField(grid)
    set!(c, 300)
    return ModifiedEVPSteppingCoefficients(c)
end

# Fallback for no specific scheme => α = β = substeps / 2
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
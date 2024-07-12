struct DynamicSteppingCoefficient{A}
    c :: A
    minimum_substeps :: Int
end

Adapt.adapt_structure(to, s::DynamicSteppingCoefficient) = 
    DynamicSteppingCoefficient(Adapt.adapt(to, s.c), s.minimum_substeps)

"""
    DynamicSteppingCoefficient(grid::AbstractGrid)

Constructs a `DynamicSteppingCoefficient` object based on the given grid that
represents the substepping coefficients for the modified EVP formulation of Kimmritz et al (2016).

In particular: α (stress stepping) == β (momentum stepping) = sqrt(4γ) 
where γ = (ζ / 4) * π² * (Δt / mᵢ) / Az is a stability parameter 

This formulation allows fast convergence in regions where sqrt(4γ) is small. Regions where
the coefficients are large correspond to regions where the ice is more solid and moves as a block
and the convergence is slower.
"""
function DynamicSteppingCoefficient(grid::AbstractGrid;
                                    minimum_substeps = 25) # smallest number of substeps
    c = CenterField(grid)
    set!(c, 300)
    return DynamicSteppingCoefficient(c, minimum_substeps)
end

@inline get_stepping_coefficients(i, j, k, grid, substeps, coefficients::DynamicSteppingCoefficient) = @inbounds coefficients.c[i, j, k] 

# Following an mEVP formulation: α = β = sqrt(4γ) (Kimmritz et al 2016)
# where γ = (ζ / 4) * π² * (Δt / mᵢ) / Az
@inline function update_stepping_coefficients!(i, j, k, grid, coefficients::DynamicSteppingCoefficient, ζ, mᵢ, Δt)
    A     = Azᶜᶜᶜ(i, j, k, grid)
    γ     = ifelse(mᵢ == 0, zero(A), ζ / 4 * π^2 * Δt / mᵢ / A)
    c_min = coefficients.minimum_substeps
    @inbounds coefficients.c[i, j, k] = max(sqrt(4γ), c_min)
    return nothing
end
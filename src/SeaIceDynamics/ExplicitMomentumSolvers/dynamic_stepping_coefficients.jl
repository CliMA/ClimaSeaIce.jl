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

In particular: α (stress stepping) == β (momentum stepping) = sqrt(γ) 
where γ = ζ * π² * (Δt / mᵢ) / Az is a stability parameter 

This formulation allows fast convergence in regions where sqrt(γ) is small. Regions where
the coefficients are large correspond to regions where the ice is more solid and moves as a block
and the convergence is slower.
"""
function DynamicSteppingCoefficient(grid::AbstractGrid;
                                    minimum_substeps = 25) # smallest number of substeps
    c = CenterField(grid)
    set!(c, 1000) # Start with a veeery large coefficient, it will be updated as soon as the viscosity is calculated

    return DynamicSteppingCoefficient(c, minimum_substeps)
end

@inline get_stepping_coefficients(i, j, k, grid, substeps, coefficients::DynamicSteppingCoefficient) = @inbounds max(coefficients.minimum_substeps, coefficients.c[i, j, k])

# Following an mEVP formulation: α = β = sqrt(γ) (Kimmritz et al 2016)
# where γ = ζ * π² * (Δt / mᵢ) / Az
@inline function update_stepping_coefficients!(i, j, k, grid, coefficients::DynamicSteppingCoefficient, ζ, mᵢ, Δt)
    A     = Azᶜᶜᶜ(i, j, k, grid)
    c_min = coefficients.minimum_substeps
    γ     = ifelse(mᵢ == 0, c_min^2, ζ * π^2 * Δt / mᵢ / A)
    @inbounds coefficients.c[i, j, k] = sqrt(γ)

    return nothing
end
struct DynamicSteppingCoefficient{A, FT}
    c :: A
    min_coeff :: FT
    max_coeff :: FT
end

Adapt.adapt_structure(to, s::DynamicSteppingCoefficient) = 
    DynamicSteppingCoefficient(Adapt.adapt(to, s.c), s.min_coeff, s.max_coeff)

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
                                    min_coeff = 50, # lower limit to `c` (as found in Kimmritz et al (2016))
                                    max_coeff = Inf)
    c = CenterField(grid)
    set!(c, max_coeff) # Start with a veeery large coefficient, it will be updated as soon as the viscosity is calculated

    FT = eltype(grid)
    return DynamicSteppingCoefficient(c, 
                                      convert(FT, min_coeff), 
                                      convert(FT, max_coeff))
end

@inline get_stepping_coefficients(i, j, k, grid, substeps, coeff::DynamicSteppingCoefficient) = @inbounds min(coeff.max_coeff, max(coeff.min_coeff, coeff.c[i, j, k]))

# Following an mEVP formulation: α = β = sqrt(γ) (Kimmritz et al 2016)
# where γ = ζ * π² * (Δt / mᵢ) / Az
@inline function update_stepping_coefficients!(i, j, k, grid, coeff::DynamicSteppingCoefficient, ζ, mᵢ, Δt)
    A     = Azᶜᶜᶜ(i, j, k, grid)
    c_max = coeff.max_coeff
    γ     = ifelse(mᵢ == 0, c_max^2, ζ * π^2 * Δt / mᵢ / A)
    @inbounds coeff.c[i, j, k] = sqrt(γ)

    return nothing
end
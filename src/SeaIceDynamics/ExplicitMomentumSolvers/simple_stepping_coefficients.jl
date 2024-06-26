
# Fallback for no specific scheme => α = β = substeps / 2
@inline get_stepping_coefficients(i, j, k, grid, substeps, ::Nothing) = substeps / 2
@inline update_stepping_coefficients!(i, j, k, grid, ::Nothing, args...) = nothing

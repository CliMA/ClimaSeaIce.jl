
# Fallback for no specific scheme => α = β = substeps 
@inline get_stepping_coefficients(i, j, k, grid, substeps, ::Nothing) = substeps 
@inline update_stepping_coefficients!(i, j, k, grid, ::Nothing, args...) = nothing

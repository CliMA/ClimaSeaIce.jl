
@inline thickness_thermodynamic_tendency(i, j, grid, h, ℵ, ::Nothing, args...) = zero(grid)

fields(::Nothing) = NamedTuple()

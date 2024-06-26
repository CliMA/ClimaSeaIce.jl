####
#### For a `Nothing` rheology we assume that u = uₒ
####

step_momentum!(model, ::Nothing, Δt, χ) = nothing

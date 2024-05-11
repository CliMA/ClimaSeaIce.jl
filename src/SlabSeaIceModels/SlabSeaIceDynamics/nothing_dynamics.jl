####
#### For a `Nothing` rheology we assume that u = uₒ
####

# The only function we need! 
step_momentum!(model, ::Nothing, Δt, χ) = nothing

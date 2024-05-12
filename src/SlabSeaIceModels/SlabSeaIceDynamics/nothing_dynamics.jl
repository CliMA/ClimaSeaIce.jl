####
#### For a `Nothing` rheology, the ice does not move!
####

# The only function we need! 
step_momentum!(model, ::Nothing, Δt, χ) = nothing

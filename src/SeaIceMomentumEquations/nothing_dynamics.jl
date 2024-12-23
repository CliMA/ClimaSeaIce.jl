####
#### For a `Nothing` dynamics, nothing happens!
####

step_momentum!(model, ::Nothing, Δt, stage) = nothing

# Fallback for nothing dynamics
compute_momentum_tendencies!(model, ::Nothing) = nothing


####
#### For a `Nothing` dynamics, nothing happens!
####

step_momentum!(model, ::Nothing, Δt, χ) = nothing

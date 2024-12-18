####
#### For `Nothing` ice dynamics, the ice velocities do not evolve 
#### -> i.e. the velocities remain equal to then initial conditions
####

# The only function we need! 
step_momentum!(model, ::Nothing, Δt) = nothing

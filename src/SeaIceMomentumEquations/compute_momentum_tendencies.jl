using Oceananigans.Utils

function compute_momentum_tendencies!(model, ::SeaIceMomentumEquations)
    
    
    u_velocity_tendency(i, j, grid,
                        clock,
                        velocities,
                        immersed_bc,
                        ocean_velocities,
                        ocean_free_surface,
                        coriolis,
                        rheology,
                        auxiliary_fields,
                        ice_thickness,
                        ice_concentration,
                        ice_density,
                        gravitational_acceleration,
                        u_top_stress,
                        u_bottom_stress,
                        u_forcing)
end


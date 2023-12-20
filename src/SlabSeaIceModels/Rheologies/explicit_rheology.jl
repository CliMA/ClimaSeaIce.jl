using Oceananigans.Grids: AbstractGrid
using Oceananigans.TimeSteppers: store_field_tendencies!
using Printf

abstract type AbstractExplicitRheology <: AbstractRheology end

function step_momentum!(model, rheology::AbstractExplicitRheology, Δt, χ)

    grid = model.grid
    arch = architecture(grid)

    u, v = model.velocities
    Guⁿ = model.timestepper.Gⁿ.u
    Gu⁻ = model.timestepper.G⁻.u
    Gvⁿ = model.timestepper.Gⁿ.v
    Gv⁻ = model.timestepper.G⁻.v

    Δτ = Δt / rheology.substeps

    for substep in 1:rheology.substeps
        
        compute_stresses!(model, model.rheology, Δτ)

        launch!(arch, grid, :xyz, 
                _compute_momentum_tendencies!, 
                model.timestepper.Gⁿ,
                grid,
                model.clock,
                model.velocities,
                model.ocean_velocities,
                model.coriolis,
                model.thickness,
                model.rheology,
                model.external_momentum_stress.u,
                model.external_momentum_stress.v,
                nothing, #model.forcing
                fields(model))

        # Advancing velocities
        launch!(arch, grid, :xyz, ab2_step_field!, u, Δτ, χ, Guⁿ, Gu⁻)
        launch!(arch, grid, :xyz, ab2_step_field!, v, Δτ, χ, Gvⁿ, Gv⁻)

        @info @sprintf("substep: %d, u: %.3e, %.3e, Gu⁻: %.3e, %.3e, Guⁿ: %.3e, %.3e", 
                       substep, extrema(u)..., extrema(Gu⁻)..., extrema(Guⁿ)...)

        # Storing tendencies
        launch!(arch, grid, :xyz, store_field_tendencies!, Gu⁻, grid, Guⁿ)
        launch!(arch, grid, :xyz, store_field_tendencies!, Gv⁻, grid, Gvⁿ)
    end

    return nothing
end

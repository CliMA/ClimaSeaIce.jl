
function step_momentum!(model, ::Nothing, Δt, χ)
    grid = model.grid
    arch = architecture(grid)

    uᵢ, vᵢ = model.velocities
    uₒ, vₒ = model.ocean_velocities

    launch!(arch, grid, :xyz, _advance_velocities!, uᵢ, vᵢ, uₒ, vₒ)

    return nothing
end

@kernel function _advance_velocities!(uᵢ, vᵢ, uₒ, vₒ)
    i, j = @index(Global, NTuple)

    @inbounds begin
        uᵢ[i, j, 1] = uₒ[i, j, 1]
        vᵢ[i, j, 1] = vₒ[i, j, 1]
    end
end

SlabSeaIceModelTendencyFields(grid, ::Nothing, tracer_names) = TracerFields(tracer_names, grid)
    
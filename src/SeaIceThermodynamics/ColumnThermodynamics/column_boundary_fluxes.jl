#####
##### Boundary contributions to the column energy solve. External heat fluxes are upward-positive: a direct flux enters
##### the top cell as -Qᵘ and the bottom cell as +Qᵇ. Dirichlet boundaries add an implicit conductance and a conductive flux.
#####

@inline function column_external_flux_temperature_slope(i, j, grid, ext, Ts, clock, model_fields)
    δT = cbrt(eps(eltype(grid))) * max(1, abs(Ts))
    Qplus = getflux(ext, i, j, grid, Ts + δT, clock, model_fields)
    Qminus = getflux(ext, i, j, grid, Ts - δT, clock, model_fields)
    return (Qplus - Qminus) / (2δT)
end

@inline function column_external_flux_temperature_slope(i, j, grid, ext::LinearizedSurfaceFlux, Ts, clock, model_fields)
    return getflux(ext.temperature_derivative, i, j, grid, Ts, clock, model_fields)
end

@inline function column_bottom_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return -Δt * λ * fields.temperature_energy_derivative[i, j, k] / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return Δt * λ * fields.temperature_energy_derivative[i, j, k] / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_flux_balance_tolerance(grid, external_flux, conductive_flux, conductance, temperature_scale)
    FT = eltype(grid)
    scale = max(1, abs(external_flux), abs(conductive_flux))
    return max(oftype(scale, 1e-6), 16eps(FT) * scale,
               8eps(FT) * abs(conductance) * max(1, abs(temperature_scale)))
end

#####
##### The surface energy balance of a MeltingConstrainedFluxBalance top
#####
##### The surface is massless: conduction from the top cell, through the ice and any snow on it, balances the
##### atmospheric flux Q(Tₛ), linearized as q + λTₛ about the last surface temperature (upward positive). A surface
##### that reaches its melting point is held there, and the heat left over melts snow or ice. Each function below
##### returns the whole surface state in one situation, and the kernel picks one:
#####
##### - `cold_column_surface`: the ice surface is below its melting point; the snow surface may be melting.
##### - `column_surface_on_held_ice`: the ice surface is held at its cap, and any snow balances on top of it.
##### - `column_surface_after_snow_melt`: the snow melts away within the step, and its latent heat joins the
#####   balance of the bare ice.
##### - `column_surface_into_reservoir`: ice thinner than the consolidation thickness passes all its surface heat to
#####   its basal reservoir.
#####
##### Each watt reaching the surface goes to exactly one of conduction Qᶜ, snow melt Qˢ and ice melt Qⁱ. C is the
##### conductance from the top cell to the ice surface, K that of the snow (zero without snow), B the two in series,
##### M the snow melt that removes all the snow in one step, and D the implicit conductance of the top row.
#####

# Temperature of a massless surface fed by conductance G from T₀ and losing q + λT + M upward. It is held at Tmax
# when conduction at Tmax brings at least the true flux there, Qmax, plus M.
@inline function held_surface_temperature(G, T₀, q, λ, M, Tmax, Qmax)
    T = min(Tmax, (G * T₀ - q - M) / (G + λ))
    return ifelse(Qmax + M <= G * (T₀ - Tmax), Tmax, T)
end

@inline function cold_column_surface(T, C, K, B, χ, q, λ, M, Tₘₛ, Qₘₛ, Q)
    Tₛ = held_surface_temperature(B, T, q, λ, zero(q), Tₘₛ, Qₘₛ)
    Tᵢ = Tₛ + (1 - 1 / (1 + χ)) * (T - Tₛ)
    Qₛ = Q(Tₛ)
    Qᶜ = C * (T - Tᵢ)
    # A melting snow surface takes the heat conduction brings, up to all the snow, and the top row sees conduction
    # alone. Otherwise the atmosphere's slope λ adds to the top row in series.
    Qˢ, D = ifelse(Tₛ >= Tₘₛ, (clamp(K * (Tᵢ - Tₛ) - Qₛ, 0, M), B), (zero(M), B * λ / (B + λ)))
    return (; Tᵢ, Tₛ, Qₛ, Qᶜ, Qˢ, Qⁱ = zero(Qᶜ), D)
end

@inline function column_surface_on_held_ice(T, C, K, R, q, λ, M, cap, Tₘₛ, Qₘₛ, Q)
    # Without snow the ice surface is the top surface.
    Tₛ = ifelse(R > 0, held_surface_temperature(K, cap, q, λ, zero(q), Tₘₛ, Qₘₛ), cap)
    Qₛ = Q(Tₛ)
    Qᶜ = C * (T - cap)
    Qˢ = ifelse(Tₛ >= Tₘₛ, clamp(K * (cap - Tₛ) - Qₛ, 0, M), zero(M))
    return (; Tᵢ = cap, Tₛ, Qₛ, Qᶜ, Qˢ, Qⁱ = max(0, Qᶜ - Qₛ - Qˢ), D = C)
end

@inline function column_surface_after_snow_melt(T, C, q, λ, M, cap, Q_cap, consolidated, Q)
    # Thin ice stays at its reservoir temperature.
    Tᵢ = ifelse(consolidated, held_surface_temperature(C, T, q, λ, M, cap, Q_cap), cap)
    Qₛ = Q(Tᵢ)
    Qᶜ = C * (T - Tᵢ)
    Qⁱ, D = ifelse(Tᵢ >= cap, (max(0, Qᶜ - Qₛ - M), C), (zero(Qᶜ), C * λ / (C + λ)))
    return (; Tᵢ, Tₛ = Tᵢ, Qₛ, Qᶜ, Qˢ = M, Qⁱ, D)
end

@inline function column_surface_into_reservoir(s)
    Qᶜ = s.Qₛ + s.Qˢ
    return (; s..., Qᶜ, Qⁱ = zero(Qᶜ), D = zero(s.D))
end

@inline column_surface_without_ice(T, s, residual) =
    (; Tᵢ = T, Tₛ = T, Qₛ = zero(s.Qₛ), Qᶜ = zero(s.Qᶜ), Qˢ = zero(s.Qˢ), Qⁱ = zero(s.Qⁱ), D = zero(s.D), residual = zero(residual))

# The surface budget and, while snow is left to conduct, the snow budget, each against its own tolerance.
@inline function column_surface_balance_residual(grid, s, T, C, K, λ, snow_conducts)
    tolerance = column_flux_balance_tolerance(grid, s.Qₛ, s.Qᶜ, C + abs(λ), max(abs(T), abs(s.Tₛ)))
    residual = abs(s.Qₛ + s.Qˢ + s.Qⁱ - s.Qᶜ) / tolerance
    snow_tolerance = column_flux_balance_tolerance(grid, s.Qₛ + s.Qˢ, K * (s.Tᵢ - s.Tₛ), K + abs(λ), max(abs(s.Tᵢ), abs(s.Tₛ)))
    snow_residual = abs(s.Qₛ + s.Qˢ - K * (s.Tᵢ - s.Tₛ)) / snow_tolerance
    return max(residual, ifelse(snow_conducts, snow_residual, zero(residual)))
end

# The snow's thermal resistance and conductance, and the melting point of the top surface: 0 °C under snow, and the
# ice surface's cap without snow, whose resistance and conductance are zero.
@inline function column_top_surface(snow, hs, h, cap)
    R = hs / (snow.conductivity * itd_factor(snow.itd_shape, h))
    return ifelse(hs > 0, (R, 1 / R, zero(cap)), (zero(hs), zero(hs), cap))
end

@kernel function _update_column_surface_fluxes!(auxiliary, fields, grid, external_flux, relation,
                                               heat_boundary_conditions, consolidation_thickness, snow, clock, model_fields, Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    @inbounds begin
        h = column_height(i, j, grid)
        T = fields.temperature[i, j, Nz]
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, Nz])
        consolidated = column_consolidated(h, consolidation_thickness[i, j, 1])
        # Unconsolidated ice exchanges heat with its isothermal basal reservoir.
        reservoir_temperature = column_bottom_reference_temperature(i, j, grid, heat_boundary_conditions.bottom, fields, relation, model_fields)
        cap = ifelse(consolidated, Tₘ, min(reservoir_temperature, Tₘ))

        hs = snow.thickness[i, j, 1]
        ρs = get_tracer(i, j, 1, grid, snow.density)
        R, K, Tₘₛ = column_top_surface(snow, hs, h, cap)
        # Backward-Euler snow balance: ρs L (hsⁿ − hs⁺) = Δt Qˢ.
        M = ρs * snow.latent_heat * hs / Δt

        C = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
        χ = C * R
        B = C / (1 + χ)

        Q(Tₛ) = getflux(external_flux, i, j, grid, Tₛ, clock, model_fields)
        Tₛ⁻ = snow.temperature[i, j, 1]
        λ = column_external_flux_temperature_slope(i, j, grid, external_flux, Tₛ⁻, clock, model_fields)
        q = Q(Tₛ⁻) - λ * Tₛ⁻
        Qₘₛ = Q(Tₘₛ)

        s = cold_column_surface(T, C, K, B, χ, q, λ, M, Tₘₛ, Qₘₛ, Q)
        on_held_ice = (s.Tᵢ >= cap) | !consolidated
        s = ifelse(on_held_ice, column_surface_on_held_ice(T, C, K, R, q, λ, M, cap, Tₘₛ, Qₘₛ, Q), s)
        snow_melts_out = (hs > 0) & (s.Tₛ >= Tₘₛ) & (K * (s.Tᵢ - s.Tₛ) - s.Qₛ >= M)
        s = ifelse(snow_melts_out, column_surface_after_snow_melt(T, C, q, λ, M, cap, Q(cap), consolidated, Q), s)
        s = ifelse(consolidated, s, column_surface_into_reservoir(s))
        residual = column_surface_balance_residual(grid, s, T, C, K, λ, (hs > 0) & !snow_melts_out)

        active = !Oceananigans.ImmersedBoundaries.immersed_cell(i, j, Nz, grid) & (h > 0)
        s = ifelse(active, (; s..., residual), column_surface_without_ice(T, s, residual))

        auxiliary.surface_temperature[i, j, 1] = s.Tᵢ
        snow.temperature[i, j, 1] = s.Tₛ
        auxiliary.surface_fluxes.heat[i, j, 1] = s.Qᶜ
        auxiliary.surface_fluxes.conductance[i, j, 1] = s.D
        auxiliary.surface_fluxes.ice_melt[i, j, 1] = s.Qⁱ
        auxiliary.surface_fluxes.snow_melt[i, j, 1] = s.Qˢ
        auxiliary.surface_fluxes.residual[i, j, 1] = s.residual
    end
end

@inline function column_bottom_boundary_energy_factor(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(i, j, k, grid, bc::MeltingConstrainedFluxBalance, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    @inbounds conductance = auxiliary.surface_fluxes.conductance[i, j, 1]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_bottom_boundary_energy_flux(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return - getflux(ext, i, j, grid, T, clock, model_fields) + λ * fields.temperature_energy_derivative[i, j, k] * fields.internal_energy[i, j, k]
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    λ = column_external_flux_temperature_slope(i, j, grid, ext, T, clock, model_fields)
    return - getflux(ext, i, j, grid, T, clock, model_fields) + λ * fields.temperature_energy_derivative[i, j, k] * fields.internal_energy[i, j, k]
end

# A Dirichlet base couples to the interior by conduction only; its external flux enters the basal Stefan balance.
@inline function column_bottom_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(i, j, grid, bc, relation, model_fields)
    return conductance * (T - Tᵇ - TE * E)
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(i, j, grid, bc, relation, model_fields)
    return conductance * (Tᵇ - T + TE * E) - getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc::MeltingConstrainedFluxBalance, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    @inbounds begin
        Q = auxiliary.surface_fluxes.heat[i, j, 1]
        C = auxiliary.surface_fluxes.conductance[i, j, 1]
        E = fields.internal_energy[i, j, k]
        TE = fields.temperature_energy_derivative[i, j, k]
    end
    return -Q + C * TE * E
end

@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    # the cap is evaluated against the start-of-step enthalpy, the same one the surface solve used
    E  = @inbounds auxiliary.surface_start_energy[i, j, 1]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    applied = min(requested, max(available, zero(available)))
    return applied - requested
end


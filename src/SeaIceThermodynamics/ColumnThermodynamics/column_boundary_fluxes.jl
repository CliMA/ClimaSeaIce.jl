#####
##### Boundary contributions to the column energy solve. External heat fluxes are upward-positive: a direct flux enters
##### the top cell as -Qᵘ and the bottom cell as +Qᵇ. Dirichlet boundaries add an implicit conductance and a conductive flux.
#####

@inline column_bottom_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)

@inline function column_external_flux_temperature_slope(i, j, grid, ext, Ts, clock, model_fields)
    δT = cbrt(eps(eltype(grid))) * max(1, abs(Ts))
    Qplus = getflux(ext, i, j, grid, Ts + δT, clock, model_fields)
    Qminus = getflux(ext, i, j, grid, Ts - δT, clock, model_fields)
    return (Qplus - Qminus) / (2δT)
end

@inline function column_external_flux_temperature_slope(i, j, grid, ext::LinearizedSurfaceFlux, Ts, clock, model_fields)
    return getflux(ext.temperature_derivative, i, j, grid, Ts, clock, model_fields)
end

@inline function column_flux_balance_tolerance(grid, external_flux, conductive_flux, conductance, temperature_scale)
    FT = eltype(grid)
    scale = max(1, abs(external_flux), abs(conductive_flux))
    return max(oftype(scale, 1e-6), 16eps(FT) * scale,
               8eps(FT) * abs(conductance) * max(1, abs(temperature_scale)))
end

@kernel function _update_column_surface_fluxes!(auxiliary, fields, grid, external_flux, relation,
                                               heat_boundary_conditions, consolidation_thickness, snow, clock, model_fields, Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)
    @inbounds begin
        T = fields.temperature[i, j, Nz]
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, fields.bulk_salinity[i, j, Nz])
        consolidated = column_consolidated(column_height(grid, i, j), consolidation_thickness[i, j, 1])
        bcs = heat_boundary_conditions
        reservoir_temperature = column_bottom_reference_temperature(bcs.bottom, i, j, grid, fields, relation)
        # Unconsolidated ice exchanges heat with its isothermal basal reservoir.
        cap_temperature = ifelse(consolidated, Tₘ, min(reservoir_temperature, Tₘ))
        Tₛ⁻ = snow.temperature[i, j, 1]
        hs = snow.thickness[i, j, 1]
        ρs = get_tracer(i, j, 1, grid, snow.density)
        R = ifelse(hs > 0, hs / (snow.conductivity * itd_factor(snow.itd_shape, column_height(grid, i, j))), zero(hs))
        # Backward-Euler snow balance: ρs L (hsⁿ − hs⁺) = Δt Qˢ.
        M = ρs * snow.latent_heat * hs / Δt
        C = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
        active = !Oceananigans.ImmersedBoundaries.immersed_cell(i, j, Nz, grid) & (column_height(grid, i, j) > 0)
        C = ifelse(active, C, one(C))
        K = ifelse(R > 0, 1 / R, zero(R))
        χ = ifelse(C == 0, zero(C), C * R)
        B = C / (1 + χ)
        Tₘₛ = ifelse(hs > 0, zero(Tₘ), cap_temperature)
        Q = getflux(external_flux, i, j, grid, Tₛ⁻, clock, model_fields)
        λ = column_external_flux_temperature_slope(i, j, grid, external_flux, Tₛ⁻, clock, model_fields)
        q = Q - λ * Tₛ⁻

        # Eliminate the massless surface from the linearized atmosphere/conduction balance.
        Tₛ = min(Tₘₛ, (B * T - q) / (B + λ))
        Qₘ = getflux(external_flux, i, j, grid, Tₘₛ, clock, model_fields)
        Tₛ = ifelse(Qₘ <= B * (T - Tₘₛ), Tₘₛ, Tₛ)
        Tᵢ = Tₛ + (1 - 1 / (1 + χ)) * (T - Tₛ)
        capped = (Tᵢ >= cap_temperature) | !consolidated
        Tᵢ = ifelse(capped, cap_temperature, Tᵢ)
        Tᶜ = ifelse(R > 0, min(Tₘₛ, (K * cap_temperature - q) / (K + λ)), cap_temperature)
        Tᶜ = ifelse(Qₘ <= K * (cap_temperature - Tₘₛ), Tₘₛ, Tᶜ)
        Tₛ = ifelse(capped, Tᶜ, Tₛ)
        Qₛ = getflux(external_flux, i, j, grid, Tₛ, clock, model_fields)
        exhausted = (hs > 0) & (Tₛ >= Tₘₛ) & (K * (Tᵢ - Tₛ) - Qₛ >= M)

        # Exhausted snow contributes its finite latent energy to the exposed-ice balance.
        Tᵦ = min(cap_temperature, (C * T - q - M) / (C + λ))
        Qᵦ = getflux(external_flux, i, j, grid, cap_temperature, clock, model_fields)
        Tᵦ = ifelse(Qᵦ + M <= C * (T - cap_temperature), cap_temperature, Tᵦ)
        Tᵢ = ifelse(exhausted & consolidated, Tᵦ, Tᵢ)
        Tₛ = ifelse(exhausted, Tᵢ, Tₛ)
        Tᵢ = ifelse(C > 0, Tᵢ, cap_temperature)
        Tₛ = ifelse((C == 0) & (hs == 0), Tᵢ, Tₛ)
        Qₛ = getflux(external_flux, i, j, grid, Tₛ, clock, model_fields)
        Qˢ = ifelse(exhausted, M, ifelse(Tₛ >= Tₘₛ, clamp(K * (Tᵢ - Tₛ) - Qₛ, 0, M), zero(M)))
        Qᶜ = ifelse(consolidated, C * (T - Tᵢ), Qₛ + Qˢ)
        Qⁱ = ifelse(consolidated & (Tᵢ >= Tₘ), max(0, Qᶜ - Qₛ - Qˢ), zero(Qᶜ))
        D = ifelse(Tₛ >= Tₘₛ, B, B * λ / (B + λ))
        D = ifelse(exhausted, C * λ / (C + λ), D)
        D = ifelse(Tᵢ >= Tₘ, C, D)
        D = ifelse(consolidated, D, zero(D))
        tolerance = column_flux_balance_tolerance(grid, Qₛ, Qᶜ, C + abs(λ), max(abs(T), abs(Tₛ)))
        residual = abs(Qₛ + Qˢ + Qⁱ - Qᶜ) / tolerance
        snow_tolerance = column_flux_balance_tolerance(grid, Qₛ + Qˢ, K * (Tᵢ - Tₛ), K + abs(λ), max(abs(Tᵢ), abs(Tₛ)))
        snow_residual = abs(Qₛ + Qˢ - K * (Tᵢ - Tₛ)) / snow_tolerance
        residual = max(residual, ifelse((hs > 0) & !exhausted, snow_residual, zero(residual)))

        auxiliary.surface_temperature[i, j, 1] = ifelse(active, Tᵢ, T)
        snow.temperature[i, j, 1] = ifelse(active, Tₛ, T)
        auxiliary.surface_fluxes.heat[i, j, 1] = ifelse(active, Qᶜ, zero(Qᶜ))
        auxiliary.surface_fluxes.conductance[i, j, 1] = ifelse(active, D, zero(D))
        auxiliary.surface_fluxes.ice_melt[i, j, 1] = ifelse(active, Qⁱ, zero(Qⁱ))
        auxiliary.surface_fluxes.snow_melt[i, j, 1] = ifelse(active, Qˢ, zero(Qˢ))
        auxiliary.surface_fluxes.residual[i, j, 1] = ifelse(active, residual, zero(residual))
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
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

# A Dirichlet base couples to the interior by conduction only; its external flux enters the basal Stefan balance.
@inline function column_bottom_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
    return conductance * (T - Tᵇ - TE * E)
end

@inline function column_top_boundary_energy_flux(i, j, k, grid, bc::ColumnDirichletBoundary, ext, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
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

# Icepack BL99 temperature-matrix solver. Used only to validate the ClimaSeaIce column thermodynamics against
# CICE/Icepack: the model itself uses the general semi-implicit energy solve, so this lives with the validation
# scripts rather than in `src`. Included by the CICE comparison driver.

using Oceananigans.Operators: Δzᶜᶜᶜ
using ClimaSeaIce.SeaIceThermodynamics:
    ColumnEnergyThermodynamics,
    FixedSalinityBrinePocketEnergyRelation,
    FluxBoundary,
    PrescribedTemperature,
    on_cpu,
    thermal_conductivity,
    ice_thermal_conductivity,
    column_dirichlet_temperature,
    getflux,
    compute_column_thermodynamic_diagnostics!,
    compute_column_shortwave_flux!,
    column_salinity_time_step!,
    internal_energy,
    melting_temperature

function solve_column_tridiagonal_vector(lower_diagonal,
                                         diagonal,
                                         upper_diagonal,
                                         rhs)
    n = length(rhs)
    solution = similar(rhs)
    scratch = zeros(eltype(rhs), max(n - 1, 0))
    β = first(diagonal)

    solution[1] = rhs[1] / β

    for k in 2:n
        scratch[k-1] = upper_diagonal[k-1] / β
        β = diagonal[k] - lower_diagonal[k-1] * scratch[k-1]
        solution[k] = (rhs[k] - lower_diagonal[k-1] * solution[k-1]) / β
    end

    for k in n-1:-1:1
        solution[k] -= scratch[k] * solution[k+1]
    end

    return solution
end

function icepack_temperature_matrix_column_step(initial_temperature,
                                                salinity,
                                                absorbed_shortwave,
                                                relation,
                                                conductivity,
                                                top_flux,
                                                bottom_temperature,
                                                Δz,
                                                Δt;
                                                max_iterations,
                                                tolerance)
    n = length(initial_temperature)
    FT = promote_type(eltype(initial_temperature), eltype(salinity), typeof(Δz), typeof(Δt))
    predicted_temperature = copy(initial_temperature)
    midpoint_conductivity =
        [ice_thermal_conductivity(conductivity, initial_temperature[k], salinity[k])
         for k in 1:n]
    interface_conductance =
        [2 * midpoint_conductivity[k-1] * midpoint_conductivity[k] /
         ((midpoint_conductivity[k-1] + midpoint_conductivity[k]) * Δz)
         for k in 2:n]
    bottom_conductance = 2 * first(midpoint_conductivity) / Δz

    phase_transitions = relation.phase_transitions
    ρᵢ = phase_transitions.density
    cᵢ = phase_transitions.heat_capacity
    L₀ = phase_transitions.reference_latent_heat
    liquidus = phase_transitions.liquidus
    top_index = n
    puny = convert(FT, 1e-11)
    previous_top_temperature_change = zero(FT)

    for iteration in 1:max_iterations
        lower_diagonal = zeros(FT, max(n - 1, 0))
        diagonal = ones(FT, n)
        upper_diagonal = zeros(FT, max(n - 1, 0))
        rhs = copy(initial_temperature)

        for k in 1:n
            Tm = melting_temperature(liquidus, salinity[k])
            heat_capacity = cᵢ - L₀ * Tm / (predicted_temperature[k] * initial_temperature[k])
            η = Δt / (ρᵢ * Δz * heat_capacity)
            rhs[k] += η * absorbed_shortwave[k]

            if n == 1
                diagonal[k] += η * bottom_conductance
                rhs[k] += η * (top_flux + bottom_conductance * bottom_temperature)
            elseif k == 1
                diagonal[k] += η * (bottom_conductance + interface_conductance[k])
                upper_diagonal[k] = -η * interface_conductance[k]
                rhs[k] += η * bottom_conductance * bottom_temperature
            elseif k == n
                lower_diagonal[k-1] = -η * interface_conductance[k-1]
                diagonal[k] += η * interface_conductance[k-1]
                rhs[k] += η * top_flux
            else
                lower_diagonal[k-1] = -η * interface_conductance[k-1]
                upper_diagonal[k] = -η * interface_conductance[k]
                diagonal[k] += η * (interface_conductance[k-1] + interface_conductance[k])
            end
        end

        next_temperature =
            solve_column_tridiagonal_vector(lower_diagonal, diagonal, upper_diagonal, rhs)

        for k in 1:n
            Tm = melting_temperature(liquidus, salinity[k])

            if Tm < 0 && next_temperature[k] > Tm - oftype(Tm, 1e-11)
                next_temperature[k] = Tm
            end
        end

        top_temperature_change = next_temperature[top_index] - predicted_temperature[top_index]
        oscillating_top_temperature =
            iteration > 1 &&
            abs(top_temperature_change) > puny &&
            abs(previous_top_temperature_change) > puny &&
            -top_temperature_change / (previous_top_temperature_change + puny^2) > 0.5

        if oscillating_top_temperature
            top_temperature_change *= 0.5

            for k in 1:n
                next_temperature[k] += 0.5 * (predicted_temperature[k] - next_temperature[k])
            end
        elseif maximum(abs.(next_temperature .- predicted_temperature)) <= tolerance
            return next_temperature
        end

        previous_top_temperature_change = top_temperature_change
        predicted_temperature = next_temperature
    end

    return predicted_temperature
end

"""
    icepack_temperature_matrix_step!(thermodynamics, external_heat_fluxes, clock, model_fields, dt;
                                     max_iterations=100, tolerance=sqrt(eps(eltype(grid))))

Advance a fixed-salinity column with the Icepack BL99 temperature-matrix
linearization. This validation utility is source-traceable to Icepack
`temperature_changes`: conductance is held at the start-of-step temperature,
and the brine-pocket heat capacity uses the secant form between the initial and
latest predicted temperatures. The supported boundary configuration is a
direct top energy flux (`FluxBoundary` with the value in `external_heat_fluxes.top`)
and a prescribed bottom temperature.
"""
function icepack_temperature_matrix_step!(thermodynamics::ColumnEnergyThermodynamics,
                                          external_heat_fluxes, clock, model_fields,
                                          Δt;
                                          max_iterations = 100,
                                          tolerance = sqrt(eps(eltype(thermodynamics.fields.internal_energy.grid))))
    grid = thermodynamics.fields.internal_energy.grid
    on_cpu(grid) || error("icepack_temperature_matrix_step! currently supports CPU grids only.")
    thermodynamics.relation isa FixedSalinityBrinePocketEnergyRelation ||
        error("icepack_temperature_matrix_step! requires FixedSalinityBrinePocketEnergyRelation.")
    thermodynamics.heat_boundary_conditions.top isa FluxBoundary ||
        error("icepack_temperature_matrix_step! requires FluxBoundary at the top boundary.")
    thermodynamics.heat_boundary_conditions.bottom isa PrescribedTemperature ||
        error("icepack_temperature_matrix_step! requires PrescribedTemperature at the bottom boundary.")

    fields = thermodynamics.fields
    auxiliary = thermodynamics.auxiliary
    relation = thermodynamics.relation
    conductivity = thermal_conductivity(thermodynamics.energy_transport)
    Nz = size(grid, 3)

    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_shortwave_flux!(thermodynamics)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            initial_temperature = [fields.temperature[i, j, k] for k in 1:Nz]
            salinity = [fields.bulk_salinity[i, j, k] for k in 1:Nz]
            absorbed_shortwave =
                [auxiliary.shortwave_flux[i, j, k+1] - auxiliary.shortwave_flux[i, j, k]
                 for k in 1:Nz]
            Δz = Δzᶜᶜᶜ(i, j, 1, grid)
            Tbot = column_dirichlet_temperature(thermodynamics.heat_boundary_conditions.bottom, i, j, grid, relation)
            # `external_heat_fluxes` are upward-positive; this BL99 reference solver adds the downward flux into the
            # top cell, so convert with a sign flip.
            top_flux = -getflux(external_heat_fluxes.top, i, j, grid, fields.temperature[i, j, Nz], clock, model_fields)

            next_temperature =
                icepack_temperature_matrix_column_step(initial_temperature,
                                                       salinity,
                                                       absorbed_shortwave,
                                                       relation,
                                                       conductivity,
                                                       top_flux,
                                                       Tbot,
                                                       Δz,
                                                       Δt;
                                                       max_iterations,
                                                       tolerance)

            for k in 1:Nz
                fields.temperature[i, j, k] = next_temperature[k]
                fields.internal_energy[i, j, k] =
                    internal_energy(relation, next_temperature[k], salinity[k])
            end
        end
    end

    compute_column_thermodynamic_diagnostics!(thermodynamics)
    column_salinity_time_step!(thermodynamics, Δt)

    return nothing
end

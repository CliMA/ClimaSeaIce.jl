using ClimaSeaIce.SeaIceThermodynamics: FluxBoundary, column_energy_time_step!
using Oceananigans
using Oceananigans.TimeSteppers: Clock
using Oceananigans.Fields: interior, ConstantField
using Oceananigans.AbstractOperations: KernelFunctionOperation
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans: fields
using Test

# Unit tests drive the column step with a model-style `external_heat_fluxes = (top, bottom)` set, a clock, no model
# fields, and a zero consolidation thickness so that every column is consolidated.
const TEST_CLOCK = Clock(time = 0.0)
const NO_FIELDS = NamedTuple()
const CONSOLIDATED = ConstantField(0.0)

insulating_fluxes() = (top = 0.0, bottom = 0.0)
insulating_boundaries() = (top = FluxBoundary(), bottom = FluxBoundary())

column_step!(thermodynamics, external_heat_fluxes, Δt; model_fields = NO_FIELDS) =
    column_energy_time_step!(thermodynamics, external_heat_fluxes, TEST_CLOCK, model_fields, CONSOLIDATED, Δt)

column_integral(c) = sum(interior(Field(c * KernelFunctionOperation{Center, Center, Center}(Δzᶜᶜᶜ, c.grid))))
column_integrated_energy(thermodynamics) = column_integral(thermodynamics.fields.internal_energy)
column_integrated_salinity(thermodynamics) = column_integral(thermodynamics.fields.bulk_salinity)

column_values(field) = vec(Array(interior(field)))

function dense_column_energy_matrix(thermodynamics)
    lower = column_values(thermodynamics.auxiliary.lower_diagonal)
    diagonal = column_values(thermodynamics.auxiliary.diagonal)
    upper = column_values(thermodynamics.auxiliary.upper_diagonal)
    Nz = length(diagonal)
    matrix = zeros(eltype(diagonal), Nz, Nz)

    for k in 1:Nz
        matrix[k, k] = diagonal[k]

        if k > 1
            matrix[k, k-1] = lower[k-1]
        end

        if k < Nz
            matrix[k, k+1] = upper[k]
        end
    end

    return matrix
end

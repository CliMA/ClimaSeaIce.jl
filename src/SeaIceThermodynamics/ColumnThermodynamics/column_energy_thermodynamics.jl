struct ColumnThermodynamicFields{E, S, T, LF, BS, TE, TS}
    internal_energy :: E
    bulk_salinity :: S
    temperature :: T
    liquid_fraction :: LF
    brine_salinity :: BS
    temperature_energy_derivative :: TE
    temperature_salinity_derivative :: TS
end

struct ColumnAuxiliaryFields{K, DE, CS, I, RHS, A, B, C, ST, BF, SE}
    thermal_conductivity :: K
    effective_energy_diffusivity :: DE
    salinity_coupling_diffusivity :: CS
    shortwave_flux :: I
    energy_rhs :: RHS
    lower_diagonal :: A
    diagonal :: B
    upper_diagonal :: C
    surface_temperature :: ST
    bottom_external_flux :: BF
    surface_start_energy :: SE
end

struct ColumnSolvers{ES, SS}
    energy_solver :: ES
    salinity_solver :: SS
end

"""
    ColumnEnergyThermodynamics(grid; kwargs...)

One-dimensional column thermodynamics with prognostic internal energy, prescribed or prognostic bulk salinity, and
diagnostic temperature, liquid fraction, brine salinity, and thermodynamic derivatives.
"""
struct ColumnEnergyThermodynamics{R, SC, ET, ST, SW, BC, F, A, SOL}
    relation :: R
    salinity_closure :: SC
    energy_transport :: ET
    salinity_transport :: ST
    shortwave_absorption :: SW
    heat_boundary_conditions :: BC
    fields :: F
    auxiliary :: A
    solvers :: SOL
end

Base.summary(::ColumnEnergyThermodynamics) = "ColumnEnergyThermodynamics"

function Base.show(io::IO, thermodynamics::ColumnEnergyThermodynamics)
    print(io, "ColumnEnergyThermodynamics", '\n')
    print(io, "├── relation: ", summary(thermodynamics.relation), '\n')
    print(io, "├── salinity_closure: ", summary(thermodynamics.salinity_closure), '\n')
    print(io, "├── energy_transport: ", summary(thermodynamics.energy_transport), '\n')
    print(io, "├── salinity_transport: ", summary(thermodynamics.salinity_transport), '\n')
    print(io, "└── shortwave_absorption: ", summary(thermodynamics.shortwave_absorption))
end

Adapt.@adapt_structure ColumnThermodynamicFields
Adapt.@adapt_structure ColumnAuxiliaryFields
Adapt.@adapt_structure ColumnEnergyThermodynamics

column_thermodynamic_fields(grid) = ColumnThermodynamicFields((Field{Center, Center, Center}(grid) for _ in 1:7)...)

function column_auxiliary_fields(grid)
    face_fields = (Field{Center, Center, Face}(grid) for _ in 1:4)
    center_fields = (Field{Center, Center, Center}(grid) for _ in 1:4)
    surface_fields = (Field{Center, Center, Nothing}(grid) for _ in 1:3)
    return ColumnAuxiliaryFields(face_fields..., center_fields..., surface_fields...)
end

settable_column_field_value(value) = value
settable_column_field_value(profile::FixedDrainedIceSalinityProfile) = z -> salinity_at_normalized_height(profile, z)

set_column_field!(field, value) = set!(field, settable_column_field_value(value))

function initialize_salinity!(fields, salinity_closure)
    profile = salinity_profile(salinity_closure)
    isnothing(profile) || set_column_field!(fields.bulk_salinity, profile)
    return nothing
end

salinity_profile(::PrognosticBulkSalinity) = nothing
salinity_profile(closure::PrescribedBulkSalinity) = closure.profile

function column_solvers(grid, auxiliary)
    solver = BatchedTridiagonalSolver(grid; auxiliary.lower_diagonal, auxiliary.diagonal, auxiliary.upper_diagonal,
                                      tridiagonal_direction = ZDirection())
    return ColumnSolvers(solver, solver)
end

default_column_heat_boundary_conditions() = (top = FluxBoundary(), bottom = FluxBoundary())

function ColumnEnergyThermodynamics(grid;
                                    relation = QuadraticLiquidusEnergyRelation(eltype(grid)),
                                    salinity_closure = PrognosticBulkSalinity(),
                                    energy_transport = ConductiveTemperatureTransport(eltype(grid)),
                                    salinity_transport = NoSalinityTransport(),
                                    shortwave_absorption = NoShortwaveAbsorption(),
                                    heat_boundary_conditions = default_column_heat_boundary_conditions(),
                                    fields = column_thermodynamic_fields(grid),
                                    auxiliary = column_auxiliary_fields(grid),
                                    solvers = column_solvers(grid, auxiliary))
    initialize_salinity!(fields, salinity_closure)
    return ColumnEnergyThermodynamics(relation, salinity_closure, energy_transport, salinity_transport, shortwave_absorption,
                                      heat_boundary_conditions, fields, auxiliary, solvers)
end

"""
    prescribed_salinity_enthalpy_thermodynamics(grid; salinity_profile = 0, kwargs...)

Construct `ColumnEnergyThermodynamics` with prescribed bulk salinity and prognostic internal energy.
"""
function prescribed_salinity_enthalpy_thermodynamics(grid; salinity_profile = 0, kw...)
    return ColumnEnergyThermodynamics(grid; salinity_closure = PrescribedBulkSalinity(salinity_profile), salinity_transport = NoSalinityTransport(), kw...)
end

"""
    evolving_salinity_mushy_thermodynamics(grid; kwargs...)

Construct `ColumnEnergyThermodynamics` with prognostic internal energy and prognostic bulk salinity.
"""
function evolving_salinity_mushy_thermodynamics(grid; energy_transport = ConductiveTemperatureTransport(eltype(grid)),
                                                salinity_transport = BulkSalinityDiffusion(eltype(grid)), kw...)
    return ColumnEnergyThermodynamics(grid; salinity_closure = PrognosticBulkSalinity(), energy_transport, salinity_transport, kw...)
end

function Oceananigans.fields(thermodynamics::ColumnEnergyThermodynamics)
    Φ = thermodynamics.fields
    return (E = Φ.internal_energy, bulk_salinity = Φ.bulk_salinity, T = Φ.temperature, liquid_fraction = Φ.liquid_fraction,
            brine_salinity = Φ.brine_salinity, temperature_energy_derivative = Φ.temperature_energy_derivative,
            temperature_salinity_derivative = Φ.temperature_salinity_derivative)
end

Oceananigans.prognostic_fields(thermodynamics::ColumnEnergyThermodynamics) = prognostic_fields(thermodynamics, thermodynamics.salinity_closure)

Oceananigans.prognostic_fields(thermodynamics::ColumnEnergyThermodynamics, ::PrescribedBulkSalinity) = (; E = thermodynamics.fields.internal_energy)

Oceananigans.prognostic_fields(thermodynamics::ColumnEnergyThermodynamics, ::PrognosticBulkSalinity) = (E = thermodynamics.fields.internal_energy, bulk_salinity = thermodynamics.fields.bulk_salinity)

function Oceananigans.Fields.set!(thermodynamics::ColumnEnergyThermodynamics; internal_energy = nothing, bulk_salinity = nothing,
                                  temperature = nothing, liquid_fraction = nothing, brine_salinity = nothing)
    Φ = thermodynamics.fields
    isnothing(internal_energy) || set_column_field!(Φ.internal_energy, internal_energy)
    isnothing(bulk_salinity)   || set_column_field!(Φ.bulk_salinity, bulk_salinity)
    isnothing(temperature)     || set_column_field!(Φ.temperature, temperature)
    isnothing(liquid_fraction) || set_column_field!(Φ.liquid_fraction, liquid_fraction)
    isnothing(brine_salinity)  || set_column_field!(Φ.brine_salinity, brine_salinity)

    energy_from_temperature = !isnothing(temperature) && isnothing(internal_energy)
    energy_from_temperature && compute_column_internal_energy!(thermodynamics)
    (energy_from_temperature || !isnothing(internal_energy)) && compute_column_thermodynamic_diagnostics!(thermodynamics)

    return nothing
end

Oceananigans.prognostic_state(thermodynamics::ColumnEnergyThermodynamics) = prognostic_state(prognostic_fields(thermodynamics))

function Oceananigans.restore_prognostic_state!(thermodynamics::ColumnEnergyThermodynamics, state)
    restore_prognostic_state!(prognostic_fields(thermodynamics), state)
    return thermodynamics
end

Oceananigans.restore_prognostic_state!(::ColumnEnergyThermodynamics, ::Nothing) = nothing


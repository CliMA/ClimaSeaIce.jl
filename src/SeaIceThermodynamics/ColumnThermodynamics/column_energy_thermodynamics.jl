struct ColumnThermodynamicFields{E, S, T, LF, BS, TE, TS}
    internal_energy :: E
    bulk_salinity :: S
    temperature :: T
    liquid_fraction :: LF
    brine_salinity :: BS
    temperature_energy_derivative :: TE
    temperature_salinity_derivative :: TS
end

struct ColumnAuxiliaryFields{K, KE, KS, DE, CS, I, RHS, A, B, C, ST, BF, SE}
    thermal_conductivity :: K
    energy_diffusivity :: KE
    salinity_diffusivity :: KS
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

One-dimensional column thermodynamics component with internal energy, bulk
salinity, temperature, liquid-fraction, brine-salinity, and thermodynamic
derivative fields. Prefer the presets
[`prescribed_salinity_enthalpy_thermodynamics`](@ref) and
[`evolving_salinity_mushy_thermodynamics`](@ref) for user-facing construction.
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
    print(io, "|-- relation: ", summary(thermodynamics.relation), '\n')
    print(io, "|-- salinity_closure: ", summary(thermodynamics.salinity_closure), '\n')
    print(io, "|-- energy_transport: ", summary(thermodynamics.energy_transport), '\n')
    print(io, "|-- salinity_transport: ", summary(thermodynamics.salinity_transport), '\n')
    print(io, "`-- shortwave_absorption: ", summary(thermodynamics.shortwave_absorption))
end

function Adapt.adapt_structure(to, fields::ColumnThermodynamicFields)
    return ColumnThermodynamicFields(Adapt.adapt(to, fields.internal_energy),
                                     Adapt.adapt(to, fields.bulk_salinity),
                                     Adapt.adapt(to, fields.temperature),
                                     Adapt.adapt(to, fields.liquid_fraction),
                                     Adapt.adapt(to, fields.brine_salinity),
                                     Adapt.adapt(to, fields.temperature_energy_derivative),
                                     Adapt.adapt(to, fields.temperature_salinity_derivative))
end

function Adapt.adapt_structure(to, auxiliary::ColumnAuxiliaryFields)
    return ColumnAuxiliaryFields(Adapt.adapt(to, auxiliary.thermal_conductivity),
                                 Adapt.adapt(to, auxiliary.energy_diffusivity),
                                 Adapt.adapt(to, auxiliary.salinity_diffusivity),
                                 Adapt.adapt(to, auxiliary.effective_energy_diffusivity),
                                 Adapt.adapt(to, auxiliary.salinity_coupling_diffusivity),
                                 Adapt.adapt(to, auxiliary.shortwave_flux),
                                 Adapt.adapt(to, auxiliary.energy_rhs),
                                 Adapt.adapt(to, auxiliary.lower_diagonal),
                                 Adapt.adapt(to, auxiliary.diagonal),
                                 Adapt.adapt(to, auxiliary.upper_diagonal),
                                 Adapt.adapt(to, auxiliary.surface_temperature),
                                 Adapt.adapt(to, auxiliary.bottom_external_flux),
                                 Adapt.adapt(to, auxiliary.surface_start_energy))
end

function Adapt.adapt_structure(to, solvers::ColumnSolvers)
    return ColumnSolvers(solvers.energy_solver, solvers.salinity_solver)
end

function Adapt.adapt_structure(to, thermodynamics::ColumnEnergyThermodynamics)
    adapted_fields = Adapt.adapt(to, thermodynamics.fields)
    adapted_auxiliary = Adapt.adapt(to, thermodynamics.auxiliary)

    return ColumnEnergyThermodynamics(Adapt.adapt(to, thermodynamics.relation),
                                      Adapt.adapt(to, thermodynamics.salinity_closure),
                                      Adapt.adapt(to, thermodynamics.energy_transport),
                                      Adapt.adapt(to, thermodynamics.salinity_transport),
                                      Adapt.adapt(to, thermodynamics.shortwave_absorption),
                                      Adapt.adapt(to, thermodynamics.heat_boundary_conditions),
                                      adapted_fields,
                                      adapted_auxiliary,
                                      Adapt.adapt(to, thermodynamics.solvers))
end

function column_thermodynamic_fields(grid)
    internal_energy = Field{Center, Center, Center}(grid)
    bulk_salinity = Field{Center, Center, Center}(grid)
    temperature = Field{Center, Center, Center}(grid)
    liquid_fraction = Field{Center, Center, Center}(grid)
    brine_salinity = Field{Center, Center, Center}(grid)
    temperature_energy_derivative = Field{Center, Center, Center}(grid)
    temperature_salinity_derivative = Field{Center, Center, Center}(grid)

    return ColumnThermodynamicFields(internal_energy,
                                     bulk_salinity,
                                     temperature,
                                     liquid_fraction,
                                     brine_salinity,
                                     temperature_energy_derivative,
                                     temperature_salinity_derivative)
end

function column_auxiliary_fields(grid)
    thermal_conductivity = Field{Center, Center, Face}(grid)
    energy_diffusivity = Field{Center, Center, Face}(grid)
    salinity_diffusivity = Field{Center, Center, Face}(grid)
    effective_energy_diffusivity = Field{Center, Center, Face}(grid)
    salinity_coupling_diffusivity = Field{Center, Center, Face}(grid)
    shortwave_flux = Field{Center, Center, Face}(grid)
    energy_rhs = Field{Center, Center, Center}(grid)
    lower_diagonal = Field{Center, Center, Center}(grid)
    diagonal = Field{Center, Center, Center}(grid)
    upper_diagonal = Field{Center, Center, Center}(grid)
    surface_temperature = Field{Center, Center, Nothing}(grid)
    bottom_external_flux = Field{Center, Center, Nothing}(grid)
    surface_start_energy = Field{Center, Center, Nothing}(grid)

    return ColumnAuxiliaryFields(thermal_conductivity,
                                 energy_diffusivity,
                                 salinity_diffusivity,
                                 effective_energy_diffusivity,
                                 salinity_coupling_diffusivity,
                                 shortwave_flux,
                                 energy_rhs,
                                 lower_diagonal,
                                 diagonal,
                                 upper_diagonal,
                                 surface_temperature,
                                 bottom_external_flux,
                                 surface_start_energy)
end

settable_column_field_value(value) = value
settable_column_field_value(profile::FixedDrainedIceSalinityProfile) =
    z -> salinity_at_normalized_height(profile, z)

set_column_field!(field, value) = set!(field, settable_column_field_value(value))

function initialize_salinity!(fields, salinity_closure)
    profile = salinity_profile(salinity_closure)
    isnothing(profile) || set_column_field!(fields.bulk_salinity, profile)
    return nothing
end

salinity_profile(::PrognosticBulkSalinity) = nothing
salinity_profile(closure::PrescribedBulkSalinity) = closure.profile

function column_solvers(grid, auxiliary)
    energy_solver = BatchedTridiagonalSolver(grid;
        lower_diagonal = auxiliary.lower_diagonal,
        diagonal = auxiliary.diagonal,
        upper_diagonal = auxiliary.upper_diagonal,
        tridiagonal_direction = ZDirection())

    return ColumnSolvers(energy_solver, energy_solver)
end

# A resting column with zero `external_heat_fluxes` is insulating on both faces.
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

    return ColumnEnergyThermodynamics(relation,
                                      salinity_closure,
                                      energy_transport,
                                      salinity_transport,
                                      shortwave_absorption,
                                      heat_boundary_conditions,
                                      fields,
                                      auxiliary,
                                      solvers)
end

"""
    prescribed_salinity_enthalpy_thermodynamics(grid; kwargs...)

Construct [`ColumnEnergyThermodynamics`](@ref) with prescribed bulk salinity
and prognostic internal energy.
"""
function prescribed_salinity_enthalpy_thermodynamics(grid;
                                                     relation = QuadraticLiquidusEnergyRelation(eltype(grid)),
                                                     salinity_profile = 0,
                                                     energy_transport = ConductiveTemperatureTransport(eltype(grid)),
                                                     shortwave_absorption = NoShortwaveAbsorption(),
                                                     heat_boundary_conditions = default_column_heat_boundary_conditions(),
                                                     kw...)
    return ColumnEnergyThermodynamics(grid;
                                      relation,
                                      salinity_closure = PrescribedBulkSalinity(salinity_profile),
                                      energy_transport,
                                      salinity_transport = NoSalinityTransport(),
                                      shortwave_absorption,
                                      heat_boundary_conditions,
                                      kw...)
end

"""
    evolving_salinity_mushy_thermodynamics(grid; kwargs...)

Construct [`ColumnEnergyThermodynamics`](@ref) with prognostic internal energy
and prognostic bulk salinity.
"""
function evolving_salinity_mushy_thermodynamics(grid;
                                                relation = QuadraticLiquidusEnergyRelation(eltype(grid)),
                                                energy_transport = ConductiveAndDiffusiveEnergyTransport(eltype(grid)),
                                                salinity_transport = BulkSalinityDiffusion(eltype(grid)),
                                                shortwave_absorption = NoShortwaveAbsorption(),
                                                heat_boundary_conditions = default_column_heat_boundary_conditions(),
                                                kw...)
    return ColumnEnergyThermodynamics(grid;
                                      relation,
                                      salinity_closure = PrognosticBulkSalinity(),
                                      energy_transport,
                                      salinity_transport,
                                      shortwave_absorption,
                                      heat_boundary_conditions,
                                      kw...)
end

function fields(thermodynamics::ColumnEnergyThermodynamics)
    thermodynamic_fields = thermodynamics.fields

    return (E = thermodynamic_fields.internal_energy,
            bulk_salinity = thermodynamic_fields.bulk_salinity,
            T = thermodynamic_fields.temperature,
            liquid_fraction = thermodynamic_fields.liquid_fraction,
            brine_salinity = thermodynamic_fields.brine_salinity,
            temperature_energy_derivative = thermodynamic_fields.temperature_energy_derivative,
            temperature_salinity_derivative = thermodynamic_fields.temperature_salinity_derivative)
end

prognostic_fields(thermodynamics::ColumnEnergyThermodynamics) = prognostic_fields(thermodynamics, thermodynamics.salinity_closure)

prognostic_fields(thermodynamics::ColumnEnergyThermodynamics, ::PrescribedBulkSalinity) =
    (; E = thermodynamics.fields.internal_energy)

prognostic_fields(thermodynamics::ColumnEnergyThermodynamics, ::PrognosticBulkSalinity) =
    (E = thermodynamics.fields.internal_energy, bulk_salinity = thermodynamics.fields.bulk_salinity)

function set!(thermodynamics::ColumnEnergyThermodynamics;
              internal_energy = nothing,
              bulk_salinity = nothing,
              temperature = nothing,
              liquid_fraction = nothing,
              brine_salinity = nothing)
    thermodynamic_fields = thermodynamics.fields

    isnothing(internal_energy) || set_column_field!(thermodynamic_fields.internal_energy, internal_energy)
    isnothing(bulk_salinity) || set_column_field!(thermodynamic_fields.bulk_salinity, bulk_salinity)
    isnothing(temperature) || set_column_field!(thermodynamic_fields.temperature, temperature)
    isnothing(liquid_fraction) || set_column_field!(thermodynamic_fields.liquid_fraction, liquid_fraction)
    isnothing(brine_salinity) || set_column_field!(thermodynamic_fields.brine_salinity, brine_salinity)

    if !isnothing(temperature) && isnothing(internal_energy)
        compute_column_internal_energy!(thermodynamics)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
    elseif !isnothing(internal_energy)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
    end

    return nothing
end

prognostic_state(thermodynamics::ColumnEnergyThermodynamics) =
    prognostic_state(prognostic_fields(thermodynamics))

function restore_prognostic_state!(thermodynamics::ColumnEnergyThermodynamics, state)
    restore_prognostic_state!(prognostic_fields(thermodynamics), state)
    return thermodynamics
end


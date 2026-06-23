import Oceananigans: prognostic_state, restore_prognostic_state!
using KernelAbstractions: @kernel, @index
using RootSolvers: SecantMethod, find_zero, CompactSolution
using Oceananigans.Architectures: CPU, architecture
using Oceananigans.Fields: AbstractField, interior
using Oceananigans.Grids: ZDirection, rnode, znode
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶜᶜᶠ, σ⁻, σⁿ
using Oceananigans.Solvers: BatchedTridiagonalSolver, solve!
using Oceananigans.Utils: launch!

@inline on_cpu(grid) = architecture(grid) isa CPU

"""
    PrescribedBulkSalinity(profile=nothing)

Closure indicating that bulk salinity is prescribed rather than prognosed.
`profile` may be any value accepted by `set!` for an Oceananigans field.
"""
struct PrescribedBulkSalinity{P}
    profile :: P
end

PrescribedBulkSalinity() = PrescribedBulkSalinity(nothing)

"""
    PrognosticBulkSalinity()

Closure indicating that bulk salinity is a prognostic column field.
"""
struct PrognosticBulkSalinity end

@inline convert_parameter(::Type{FT}, value::Number) where FT = convert(FT, value)
@inline convert_parameter(::Type{FT}, value) where FT = value

"""
    MaykutUntersteinerConductivity([FT=Oceananigans.defaults.FloatType; kwargs...])

CICE/Icepack BL99 `conduct = "MU71"` thermal conductivity closure.
"""
struct MaykutUntersteinerConductivity{FT}
    fresh_ice_conductivity :: FT
    salinity_coefficient :: FT
    minimum_conductivity :: FT
    temperature_floor :: FT
end

MaykutUntersteinerConductivity(FT::DataType=Oceananigans.defaults.FloatType;
                               fresh_ice_conductivity = 2.03,
                               salinity_coefficient = 0.13,
                               minimum_conductivity = 0.10,
                               temperature_floor = 1e-11) =
    MaykutUntersteinerConductivity(convert(FT, fresh_ice_conductivity),
                                   convert(FT, salinity_coefficient),
                                   convert(FT, minimum_conductivity),
                                   convert(FT, temperature_floor))

Base.summary(::MaykutUntersteinerConductivity) = "MaykutUntersteinerConductivity"

function Adapt.adapt_structure(to, conductivity::MaykutUntersteinerConductivity)
    return MaykutUntersteinerConductivity(Adapt.adapt(to, conductivity.fresh_ice_conductivity),
                                          Adapt.adapt(to, conductivity.salinity_coefficient),
                                          Adapt.adapt(to, conductivity.minimum_conductivity),
                                          Adapt.adapt(to, conductivity.temperature_floor))
end

"""
    BubblyBrineConductivity([FT=Oceananigans.defaults.FloatType; kwargs...])

CICE/Icepack BL99 `conduct = "bubbly"` thermal conductivity closure.
"""
struct BubblyBrineConductivity{FT}
    ice_density :: FT
    pure_ice_density :: FT
    minimum_conductivity :: FT
    temperature_floor :: FT
end

BubblyBrineConductivity(FT::DataType=Oceananigans.defaults.FloatType;
                        ice_density = 917,
                        pure_ice_density = 917,
                        minimum_conductivity = 0.10,
                        temperature_floor = 1e-11) =
    BubblyBrineConductivity(convert(FT, ice_density),
                            convert(FT, pure_ice_density),
                            convert(FT, minimum_conductivity),
                            convert(FT, temperature_floor))

Base.summary(::BubblyBrineConductivity) = "BubblyBrineConductivity"

function Adapt.adapt_structure(to, conductivity::BubblyBrineConductivity)
    return BubblyBrineConductivity(Adapt.adapt(to, conductivity.ice_density),
                                   Adapt.adapt(to, conductivity.pure_ice_density),
                                   Adapt.adapt(to, conductivity.minimum_conductivity),
                                   Adapt.adapt(to, conductivity.temperature_floor))
end

@inline ice_thermal_conductivity(conductivity::Number, temperature, bulk_salinity) = conductivity

@inline function ice_thermal_conductivity(conductivity::MaykutUntersteinerConductivity,
                                          temperature,
                                          bulk_salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.fresh_ice_conductivity +
        conductivity.salinity_coefficient * bulk_salinity / T
    return max(K, conductivity.minimum_conductivity)
end

@inline function ice_thermal_conductivity(conductivity::BubblyBrineConductivity,
                                          temperature,
                                          bulk_salinity)
    T = min(-conductivity.temperature_floor, temperature)
    K = conductivity.ice_density / conductivity.pure_ice_density *
        (2.11 - 0.011 * temperature + 0.09 * bulk_salinity / T)
    return max(K, conductivity.minimum_conductivity)
end

@inline function face_thermal_conductivity(conductivity,
                                           lower_temperature,
                                           lower_salinity,
                                           lower_thickness,
                                           upper_temperature,
                                           upper_salinity,
                                           upper_thickness)
    lower_conductivity = ice_thermal_conductivity(conductivity,
                                                  lower_temperature,
                                                  lower_salinity)

    upper_conductivity = ice_thermal_conductivity(conductivity,
                                                  upper_temperature,
                                                  upper_salinity)

    return (lower_thickness + upper_thickness) /
           (lower_thickness / lower_conductivity +
            upper_thickness / upper_conductivity)
end

"""
    ConductiveTemperatureTransport([FT=Oceananigans.defaults.FloatType; conductivity=2])

Energy transport closure for conductive heat flux proportional to the
temperature gradient.
"""
struct ConductiveTemperatureTransport{K}
    conductivity :: K
end

ConductiveTemperatureTransport(FT::DataType=Oceananigans.defaults.FloatType;
                               conductivity = 2) =
    ConductiveTemperatureTransport(convert_parameter(FT, conductivity))

"""
    DiffusiveEnergyTransport([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Energy transport closure for direct diffusion of internal energy.
"""
struct DiffusiveEnergyTransport{K}
    diffusivity :: K
end

DiffusiveEnergyTransport(FT::DataType=Oceananigans.defaults.FloatType;
                         diffusivity = 0) =
    DiffusiveEnergyTransport(convert(FT, diffusivity))

"""
    ConductiveAndDiffusiveEnergyTransport([FT=Oceananigans.defaults.FloatType; conductivity=2, diffusivity=0])

Energy transport closure combining conductive temperature transport and direct
internal-energy diffusion.
"""
struct ConductiveAndDiffusiveEnergyTransport{K, D}
    conductivity :: K
    diffusivity :: D
end

ConductiveAndDiffusiveEnergyTransport(FT::DataType=Oceananigans.defaults.FloatType;
                                      conductivity = 2,
                                      diffusivity = 0) =
    ConductiveAndDiffusiveEnergyTransport(convert_parameter(FT, conductivity),
                                          convert(FT, diffusivity))

"""
    NoSalinityTransport()

Salinity transport closure that leaves prognostic bulk salinity unchanged.
"""
struct NoSalinityTransport end

"""
    BulkSalinityDiffusion([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Closed-boundary scalar diffusion closure for prognostic bulk salinity.
"""
struct BulkSalinityDiffusion{K}
    diffusivity :: K
end

BulkSalinityDiffusion(FT::DataType=Oceananigans.defaults.FloatType;
                      diffusivity = 0) =
    BulkSalinityDiffusion(convert(FT, diffusivity))

"""
    BrineSalinityDiffusion([FT=Oceananigans.defaults.FloatType; diffusivity=0])

Marker closure for brine-salinity diffusion. The first column implementation
uses [`BulkSalinityDiffusion`](@ref) for the scalar salinity step.
"""
struct BrineSalinityDiffusion{K}
    diffusivity :: K
end

BrineSalinityDiffusion(FT::DataType=Oceananigans.defaults.FloatType;
                       diffusivity = 0) =
    BrineSalinityDiffusion(convert(FT, diffusivity))

"""
    NoShortwaveAbsorption()

Shortwave-radiation closure with no internal shortwave source.
"""
struct NoShortwaveAbsorption end

"""
    ExponentialShortwaveAbsorption([FT=Oceananigans.defaults.FloatType; surface_transmission=0, attenuation_scale=1])

Beer-law shortwave-radiation closure. `surface_transmission` is the transmitted
shortwave flux at the top face of the column and `attenuation_scale` is the
vertical e-folding length.
"""
struct ExponentialShortwaveAbsorption{T, L}
    surface_transmission :: T
    attenuation_scale :: L
end

ExponentialShortwaveAbsorption(FT::DataType=Oceananigans.defaults.FloatType;
                               surface_transmission = 0,
                               attenuation_scale = 1) =
    ExponentialShortwaveAbsorption(convert(FT, surface_transmission),
                                   convert(FT, attenuation_scale))

struct ColumnThermodynamicFields{E, S, T, LF, BS, TE, TS}
    internal_energy :: E
    bulk_salinity :: S
    temperature :: T
    liquid_fraction :: LF
    brine_salinity :: BS
    temperature_energy_derivative :: TE
    temperature_salinity_derivative :: TS
end

struct ColumnAuxiliaryFields{K, KE, KS, DE, CS, I, RHS, A, B, C, ST}
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
                                 Adapt.adapt(to, auxiliary.surface_temperature))
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
                                 surface_temperature)
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

prognostic_fields(thermodynamics::ColumnEnergyThermodynamics) =
    prognostic_fields(thermodynamics, thermodynamics.salinity_closure)

function prognostic_fields(thermodynamics::ColumnEnergyThermodynamics,
                           ::PrescribedBulkSalinity)
    return (; E = thermodynamics.fields.internal_energy)
end

function prognostic_fields(thermodynamics::ColumnEnergyThermodynamics,
                           ::PrognosticBulkSalinity)
    return (E = thermodynamics.fields.internal_energy,
            bulk_salinity = thermodynamics.fields.bulk_salinity)
end

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

#####
##### Diagnostics and column energy step
#####

thermal_conductivity(transport::ConductiveTemperatureTransport) = transport.conductivity
thermal_conductivity(::DiffusiveEnergyTransport) = 0
thermal_conductivity(transport::ConductiveAndDiffusiveEnergyTransport) = transport.conductivity

energy_diffusivity(::ConductiveTemperatureTransport) = 0
energy_diffusivity(transport::DiffusiveEnergyTransport) = transport.diffusivity
energy_diffusivity(transport::ConductiveAndDiffusiveEnergyTransport) = transport.diffusivity

salinity_diffusivity(::NoSalinityTransport) = 0
salinity_diffusivity(transport::BulkSalinityDiffusion) = transport.diffusivity

@inline shortwave_flux(::NoShortwaveAbsorption, z, z_top) = zero(z)

@inline function shortwave_flux(absorption::ExponentialShortwaveAbsorption, z, z_top)
    return absorption.surface_transmission *
           exp((z - z_top) / absorption.attenuation_scale)
end

@kernel function _compute_column_internal_energy!(fields, relation)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        fields.internal_energy[i, j, k] =
            internal_energy(relation,
                            fields.temperature[i, j, k],
                            fields.bulk_salinity[i, j, k])
    end
end

"""
    compute_column_internal_energy!(thermodynamics)

Fill the internal-energy field from the current temperature and bulk-salinity
fields.
"""
function compute_column_internal_energy!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_internal_energy!,
            thermodynamics.fields,
            thermodynamics.relation)

    return nothing
end

@kernel function _compute_column_thermodynamic_diagnostics!(fields, relation)
    i, j, k = @index(Global, NTuple)

    @inbounds begin
        E = fields.internal_energy[i, j, k]
        S = fields.bulk_salinity[i, j, k]

        fields.temperature[i, j, k] =
            temperature(relation, E, S)

        fields.liquid_fraction[i, j, k] =
            liquid_fraction(relation, E, S)

        fields.brine_salinity[i, j, k] =
            brine_salinity(relation, E, S)

        fields.temperature_energy_derivative[i, j, k] =
            temperature_energy_derivative(relation, E, S)

        fields.temperature_salinity_derivative[i, j, k] =
            temperature_salinity_derivative(relation, E, S)
    end
end

"""
    compute_column_thermodynamic_diagnostics!(thermodynamics)

Fill temperature, liquid fraction, brine salinity, and thermodynamic derivative
diagnostics from internal energy and bulk salinity.
"""
function compute_column_thermodynamic_diagnostics!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_thermodynamic_diagnostics!,
            thermodynamics.fields,
            thermodynamics.relation)

    return nothing
end

@kernel function _compute_column_transport_coefficients!(auxiliary, fields, grid, energy_transport)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    kappa_T = thermal_conductivity(energy_transport)
    kappa_E = energy_diffusivity(energy_transport)

    @inbounds begin
        T = fields.temperature
        S = fields.bulk_salinity

        if k == 1
            boundary_conductivity = ice_thermal_conductivity(kappa_T, T[i, j, 1], S[i, j, 1])
            TE = fields.temperature_energy_derivative[i, j, 1]
            TS = fields.temperature_salinity_derivative[i, j, 1]

            auxiliary.thermal_conductivity[i, j, 1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, 1] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, 1] = boundary_conductivity * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, 1] = boundary_conductivity * TS
        end

        if k > 1
            TE = (fields.temperature_energy_derivative[i, j, k-1] +
                  fields.temperature_energy_derivative[i, j, k]) / 2

            TS = (fields.temperature_salinity_derivative[i, j, k-1] +
                  fields.temperature_salinity_derivative[i, j, k]) / 2

            face_conductivity =
                face_thermal_conductivity(kappa_T,
                                          T[i, j, k-1],
                                          S[i, j, k-1],
                                          Δzᶜᶜᶜ(i, j, k-1, grid),
                                          T[i, j, k],
                                          S[i, j, k],
                                          Δzᶜᶜᶜ(i, j, k, grid))

            auxiliary.thermal_conductivity[i, j, k] = face_conductivity
            auxiliary.energy_diffusivity[i, j, k] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, k] = face_conductivity * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, k] = face_conductivity * TS
        end

        if k == Nz
            kp1 = Nz + 1
            boundary_conductivity = ice_thermal_conductivity(kappa_T, T[i, j, Nz], S[i, j, Nz])
            TE = fields.temperature_energy_derivative[i, j, Nz]
            TS = fields.temperature_salinity_derivative[i, j, Nz]

            auxiliary.thermal_conductivity[i, j, kp1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, kp1] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, kp1] = boundary_conductivity * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, kp1] = boundary_conductivity * TS
        end
    end
end

"""
    compute_column_transport_coefficients!(thermodynamics)

Compute face-centered energy transport coefficients for the semi-implicit
energy solve.
"""
function compute_column_transport_coefficients!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_transport_coefficients!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.energy_transport)

    return nothing
end

@kernel function _compute_column_salinity_diffusivity!(auxiliary, fields, salinity_transport)
    i, j, k = @index(Global, NTuple)
    Nz = size(fields.bulk_salinity.grid, 3)

    kappa_S = salinity_diffusivity(salinity_transport)

    @inbounds begin
        if k == 1
            auxiliary.salinity_diffusivity[i, j, 1] = zero(kappa_S)
        end

        if k > 1
            auxiliary.salinity_diffusivity[i, j, k] = kappa_S
        end

        if k == Nz
            auxiliary.salinity_diffusivity[i, j, Nz+1] = zero(kappa_S)
        end
    end
end

"""
    compute_column_salinity_diffusivity!(thermodynamics)

Compute face-centered bulk-salinity diffusivity for the scalar salinity solve.
"""
function compute_column_salinity_diffusivity!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_salinity_diffusivity!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            thermodynamics.salinity_transport)

    return nothing
end

@kernel function _compute_column_shortwave_flux!(auxiliary, grid, shortwave_absorption)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)
    z_top = znode(i, j, Nz+1, grid, Center(), Center(), Face())

    @inbounds begin
        if k == 1
            z_bottom = znode(i, j, 1, grid, Center(), Center(), Face())
            auxiliary.shortwave_flux[i, j, 1] =
                shortwave_flux(shortwave_absorption, z_bottom, z_top)
        end

        if k > 1
            z_face = znode(i, j, k, grid, Center(), Center(), Face())
            auxiliary.shortwave_flux[i, j, k] =
                shortwave_flux(shortwave_absorption, z_face, z_top)
        end

        if k == Nz
            auxiliary.shortwave_flux[i, j, Nz+1] =
                shortwave_flux(shortwave_absorption, z_top, z_top)
        end
    end
end

"""
    compute_column_shortwave_flux!(thermodynamics)

Compute face-centered shortwave fluxes for the current shortwave absorption
closure.
"""
function compute_column_shortwave_flux!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _compute_column_shortwave_flux!,
            thermodynamics.auxiliary,
            grid,
            thermodynamics.shortwave_absorption)

    return nothing
end

@kernel function _compute_column_surface_stefan_residual_flux!(residual_flux,
                                                               auxiliary,
                                                               fields,
                                                               grid,
                                                               heat_boundary_conditions,
                                                               external_heat_fluxes,
                                                               relation,
                                                               clock,
                                                               model_fields,
                                                               Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds residual_flux[i, j, 1] =
        column_surface_stefan_residual_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                            i, j, Nz, grid, auxiliary, fields, relation, clock, model_fields, Δt)
end

"""
    compute_column_surface_stefan_residual_flux!(residual_flux, thermodynamics, external_heat_fluxes, clock, model_fields, dt)

Fill `residual_flux` with the top-surface Stefan residual implied by a
`MeltingConstrainedFluxBalance` surface boundary. The residual is positive for
growth and negative for melt, matching [`column_stefan_thickness_update!`](@ref).
"""
function compute_column_surface_stefan_residual_flux!(residual_flux,
                                                      thermodynamics::ColumnEnergyThermodynamics,
                                                      external_heat_fluxes,
                                                      clock,
                                                      model_fields,
                                                      Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _compute_column_surface_stefan_residual_flux!,
            residual_flux,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.heat_boundary_conditions,
            external_heat_fluxes,
            thermodynamics.relation,
            clock,
            model_fields,
            Δt)

    return nothing
end

@inline diffusion_factor(i, j, kᶜ, kᶠ, grid, K, Δt) =  Δt * @inbounds(K[i, j, kᶠ]) / (Δzᶜᶜᶜ(i, j, kᶜ, grid) * Δzᶜᶜᶠ(i, j, kᶠ, grid))

# Conservative metric Jacobian: the column-thickness ratio (σ ≡ column height / reference height).
@inline previous_to_current_column_metric_ratio(i, j, k, grid::SeaIceColumnGrid) =
    previous_column_height(grid, i, j) / column_height(grid, i, j)

# The two interfaces move independently, so the swept-face displacement is the base motion plus the standard
# thickness-stretch term (σⁿ - σ⁻) r: basal growth/melt moves only the base, surface melt/snow-ice only the surface.
@inline function column_face_displacement(i, j, k, grid::SeaIceColumnGrid)
    z = sea_ice_discretization(grid)
    δbase = @inbounds z.hbⁿ[i, j, 1] - z.hb⁻[i, j, 1]
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Face())
    σᶜᶜ⁻ = σ⁻(i, j, k, grid, Center(), Center(), Face())
    r = rnode(i, j, k, grid, Center(), Center(), Face())
    return δbase + (σᶜᶜⁿ - σᶜᶜ⁻) * r
end

# The temperature of a Dirichlet boundary (ocean equilibrium or prescribed), shared with the slab vocabulary.
@inline column_dirichlet_temperature(bc, i, j, grid, relation) =
    bottom_temperature(i, j, grid, bc, relation.phase_transitions.liquidus)

# Swept material on a moving interface takes the adjacent interior enthalpy by default. A Dirichlet basal interface
# (congelation) instead lays down ice at the interface temperature — the CICE BL99 new-ice enthalpy.
@inline moving_bottom_boundary_value(boundary, i, j, grid, fields, relation, default) = default
@inline moving_top_boundary_value(boundary, i, j, grid, fields, relation, default) = default

@inline function moving_bottom_boundary_value(boundary::Union{PrescribedTemperature, IceWaterThermalEquilibrium}, i, j, grid, fields, relation, default)
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    S = @inbounds fields.bulk_salinity[i, j, 1]
    return internal_energy(relation, Tᵇ, S)
end

@inline function moving_face_value(i, j, k, grid, field, face_displacement, bottom_boundary_value, top_boundary_value)
    Nz = size(grid, 3)

    if k == 1
        interior_value = @inbounds field[i, j, 1]
        return face_displacement < 0 ? bottom_boundary_value : interior_value
    elseif k == Nz + 1
        interior_value = @inbounds field[i, j, Nz]
        return face_displacement > 0 ? top_boundary_value : interior_value
    elseif face_displacement >= 0
        return @inbounds field[i, j, k]
    else
        return @inbounds field[i, j, k-1]
    end
end

@inline function moving_face_displacement_flux(i, j, k, grid, field)
    displacement = column_face_displacement(i, j, k, grid)
    Nz = size(grid, 3)
    bottom_boundary_value = @inbounds field[i, j, 1]
    top_boundary_value = @inbounds field[i, j, Nz]
    value = moving_face_value(i, j, k, grid, field, displacement,
                              bottom_boundary_value, top_boundary_value)
    return displacement * value
end

@inline function moving_face_displacement_flux(i, j, k, grid, field, fields, relation,
                                               bottom_boundary,
                                               top_boundary)
    displacement = column_face_displacement(i, j, k, grid)
    Nz = size(grid, 3)
    bottom_default = @inbounds field[i, j, 1]
    top_default = @inbounds field[i, j, Nz]
    bottom_boundary_value = moving_bottom_boundary_value(bottom_boundary, i, j, grid, fields, relation, bottom_default)
    top_boundary_value = moving_top_boundary_value(top_boundary, i, j, grid, fields, relation, top_default)
    value = moving_face_value(i, j, k, grid, field, displacement, bottom_boundary_value, top_boundary_value)
    return displacement * value
end

@inline salinity_coupling_flux(i, j, k, grid, C, S) = @inbounds(C[i, j, k] * (S[i, j, k] - S[i, j, k-1]) / Δzᶜᶜᶠ(i, j, k, grid))

@inline column_boundary_temperature_conductance(i, j, kf, kc, grid, auxiliary) = @inbounds 2 * auxiliary.thermal_conductivity[i, j, kf] / Δzᶜᶜᶜ(i, j, kc, grid)

@inline column_consolidated(consolidation_thickness::Nothing, grid, i, j) = true
@inline function column_consolidated(consolidation_thickness, grid, i, j)
    h = column_height(grid, i, j)
    hᶜ = @inbounds consolidation_thickness[i, j, 1]
    return h ≥ hᶜ
end

@inline function column_energy_to_complete_melt(relation, E, S, Δz)
    return (complete_melt_energy(relation, S) - E) * Δz
end

#####
##### Boundary contributions to the column energy solve, keyed on (heat_boundary_condition, external_flux).
##### `external_heat_fluxes` use the same upward-positive convention as `SlabThermodynamics` (positive = heat
##### leaving the ice, i.e. into the air at the top and from the ocean into the base at the bottom), so the same
##### `model.external_heat_fluxes` drives a slab and a column identically. A direct flux therefore enters the top
##### cell as `-Qᵘ` and the bottom cell as `+Qᵇ`. Dirichlet boundaries (`PrescribedTemperature`,
##### `IceWaterThermalEquilibrium`) additionally add an implicit conductance factor (LHS) and a conductive flux.
#####

const ColumnDirichletBoundary = Union{PrescribedTemperature, IceWaterThermalEquilibrium}

@inline column_bottom_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(Δt)

@inline function column_bottom_boundary_energy_factor(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / Δzᶜᶜᶜ(i, j, k, grid)
end

# A direct-flux boundary injects its paired upward-positive external flux: a `left_flux` of `-Qᵇ` adds ocean heat
# to the base, a `right_flux` of `-Qᵘ` removes heat to the air at the top.
@inline function column_bottom_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    T = @inbounds fields.temperature[i, j, k]
    return - getflux(ext, i, j, grid, T, clock, model_fields)
end

@inline function column_bottom_boundary_energy_flux(bc::ColumnDirichletBoundary,
                                                    ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
    return conductance * (T - Tᵇ - TE * E) - getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
end

@inline function column_top_boundary_energy_flux(bc::ColumnDirichletBoundary, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T  = @inbounds fields.temperature[i, j, k]
    E  = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    Tᵇ = column_dirichlet_temperature(bc, i, j, grid, relation)
    return conductance * (Tᵇ - T + TE * E) - getflux(ext, i, j, grid, Tᵇ, clock, model_fields)
end

# A `MeltingConstrainedFluxBalance` surface injects the downward flux `-Qᵘ(Tₛ)` evaluated at the iteratively-solved
# surface temperature `Tₛ`, capped at the top-cell complete-melt energy; the excess becomes the surface Stefan
# residual that drives surface melt.
@inline function column_top_boundary_energy_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    E  = @inbounds fields.internal_energy[i, j, k]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    return min(requested, max(available, zero(available)))
end

# The surface Stefan residual is the part of the requested top flux not absorbed by the capped cell warming.
@inline column_surface_stefan_residual_flux(bc, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt) = zero(eltype(grid))

@inline function column_surface_stefan_residual_flux(bc::MeltingConstrainedFluxBalance, ext, i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
    Δz = Δzᶜᶜᶜ(i, j, k, grid)
    E  = @inbounds fields.internal_energy[i, j, k]
    S  = @inbounds fields.bulk_salinity[i, j, k]
    Tₛ = @inbounds auxiliary.surface_temperature[i, j, 1]
    requested = - getflux(ext, i, j, grid, Tₛ, clock, model_fields)
    available = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    applied = min(requested, max(available, zero(available)))
    return applied - requested
end

@kernel function _assemble_column_energy_system!(auxiliary, fields, grid, heat_boundary_conditions, external_heat_fluxes, relation, consolidation_thickness, clock, model_fields, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    if column_consolidated(consolidation_thickness, grid, i, j)
    @inbounds begin
        D = auxiliary.effective_energy_diffusivity
        C = auxiliary.salinity_coupling_diffusivity
        I = auxiliary.shortwave_flux
        S = fields.bulk_salinity

        left_factor = if k == 1
            zero(Δt)
        else
            diffusion_factor(i, j, k, k, grid, D, Δt)
        end

        right_factor = if k == Nz
            zero(Δt)
        else
            diffusion_factor(i, j, k, k+1, grid, D, Δt)
        end

        bottom_factor = if k == 1
            column_bottom_boundary_energy_factor(heat_boundary_conditions.bottom, external_heat_fluxes.bottom,
                                                 i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            zero(Δt)
        end

        top_factor = if k == Nz
            column_top_boundary_energy_factor(heat_boundary_conditions.top, external_heat_fluxes.top,
                                              i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            zero(Δt)
        end

        auxiliary.diagonal[i, j, k] = 1 + left_factor + right_factor + bottom_factor + top_factor

        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = -left_factor
        end

        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = -right_factor
        else
            auxiliary.upper_diagonal[i, j, k] = zero(right_factor)
            auxiliary.lower_diagonal[i, j, k] = zero(left_factor)
        end

        left_flux = if k == 1
            column_bottom_boundary_energy_flux(heat_boundary_conditions.bottom, external_heat_fluxes.bottom,
                                               i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            salinity_coupling_flux(i, j, k, grid, C, S)
        end
        left_flux += I[i, j, k]

        right_flux = if k == Nz
            column_top_boundary_energy_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                            i, j, k, grid, auxiliary, fields, relation, clock, model_fields, Δt)
        else
            salinity_coupling_flux(i, j, k+1, grid, C, S)
        end
        right_flux += I[i, j, k+1]

        flux_tendency = (right_flux - left_flux) / Δzᶜᶜᶜ(i, j, k, grid)

        # Conservative moving-grid balance divided by the current layer thickness:
        # Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ).
        moving_flux_difference =
            moving_face_displacement_flux(i, j, k+1, grid, fields.internal_energy, fields, relation,
                                          heat_boundary_conditions.bottom,
                                          heat_boundary_conditions.top) -
            moving_face_displacement_flux(i, j, k, grid, fields.internal_energy, fields, relation,
                                          heat_boundary_conditions.bottom,
                                          heat_boundary_conditions.top)
        moving_tendency = moving_flux_difference / Δzᶜᶜᶜ(i, j, k, grid)
        metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] =
            metric_ratio * fields.internal_energy[i, j, k] + moving_tendency + Δt * flux_tendency
    end
    else
    @inbounds begin
        # Unconsolidated column: no resolved interior profile. Hold the internal energy with an identity row;
        # the thickness change is taken up by the slab flux balance in the volume update.
        auxiliary.diagonal[i, j, k] = one(Δt)
        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = zero(Δt)
        end
        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = zero(Δt)
        else
            auxiliary.lower_diagonal[i, j, k] = zero(Δt)
        end
        auxiliary.energy_rhs[i, j, k] = fields.internal_energy[i, j, k]
    end
    end
end

"""
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, dt)

Assemble the tridiagonal backward-Euler system for one internal-energy step.
`external_heat_fluxes` is the model's `(top, bottom)` flux set, evaluated through
`getflux` with `clock` and `model_fields`. A `consolidation_thickness` of `nothing`
treats every column as consolidated (the standalone path); otherwise columns thinner
than it get an identity row (no resolved profile).
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics,
                                        external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    compute_column_shortwave_flux!(thermodynamics)

    launch!(arch, grid, :xyz,
            _assemble_column_energy_system!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            thermodynamics.heat_boundary_conditions,
            external_heat_fluxes,
            thermodynamics.relation,
            consolidation_thickness,
            clock,
            model_fields,
            Δt)

    return nothing
end

#####
##### Implicit surface-temperature solve for a `MeltingConstrainedFluxBalance` top, Icepack `temperature_changes`
##### style: an outer iteration around the tridiagonal solve makes the lagged top flux self-consistent with the
##### surface temperature `Tₛ`, found from the massless-surface balance `Q_atm(Tₛ) + conductance·(T_top − Tₛ) = 0`
##### and capped at melting. Negligible conduction decouples the surface, so `Tₛ` is left unchanged there.
#####

@inline function column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
    FT = eltype(grid)
    conductance > eps(FT) || return Tₛ⁻
    # Massless-surface balance in the upward-positive convention: Qᵘ(Tₛ) = conductance·(T_top − Tₛ).
    balance(T) = getflux(external_flux, i, j, grid, T, clock, model_fields) - conductance * (T_top - T)
    solution = find_zero(balance, SecantMethod{FT}(Tₛ⁻ + one(FT), Tₛ⁻), CompactSolution())
    return solution.root
end

@kernel function _update_column_surface_temperature!(auxiliary, fields, grid, external_flux, relation, clock, model_fields)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        Tₛ⁻ = auxiliary.surface_temperature[i, j, 1]
        T_top = fields.temperature[i, j, Nz]
        S_top = fields.bulk_salinity[i, j, Nz]
        conductance = column_boundary_temperature_conductance(i, j, Nz+1, Nz, grid, auxiliary)
        Tₘ = melting_temperature(relation.phase_transitions.liquidus, S_top)
        Tₛ = column_surface_temperature_balance(i, j, grid, Tₛ⁻, conductance, T_top, external_flux, clock, model_fields)
        auxiliary.surface_temperature[i, j, 1] = min(Tₛ, Tₘ)
    end
end

"""
    solve_column_energy_system!(thermodynamics)

Solve the assembled column energy system into the internal-energy field.
"""
function solve_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics)
    solve!(thermodynamics.fields.internal_energy,
           thermodynamics.solvers.energy_solver,
           thermodynamics.auxiliary.energy_rhs)

    return nothing
end

@kernel function _assemble_column_salinity_system!(auxiliary, fields, grid, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        K = auxiliary.salinity_diffusivity

        left_factor = if k == 1
            zero(Δt)
        else
            diffusion_factor(i, j, k, k, grid, K, Δt)
        end

        right_factor = if k == Nz
            zero(Δt)
        else
            diffusion_factor(i, j, k, k+1, grid, K, Δt)
        end

        auxiliary.diagonal[i, j, k] = 1 + left_factor + right_factor

        if k > 1
            auxiliary.lower_diagonal[i, j, k-1] = -left_factor
        end

        if k < Nz
            auxiliary.upper_diagonal[i, j, k] = -right_factor
        else
            auxiliary.upper_diagonal[i, j, k] = zero(right_factor)
            auxiliary.lower_diagonal[i, j, k] = zero(left_factor)
        end

        metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)

        # Same conservative moving-grid update as energy, without boundary heat
        # fluxes: Jⁿ⁺¹ Δr Sⁿ⁺¹ = Jⁿ Δr Sⁿ + δzᵤ Sᵘᵖᵤ - δzₗ Sᵘᵖₗ.
        moving_flux_difference =
            moving_face_displacement_flux(i, j, k+1, grid, fields.bulk_salinity) -
            moving_face_displacement_flux(i, j, k, grid, fields.bulk_salinity)
        moving_tendency = moving_flux_difference / Δzᶜᶜᶜ(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] = metric_ratio * fields.bulk_salinity[i, j, k] + moving_tendency
    end
end

"""
    assemble_column_salinity_system!(thermodynamics, dt)

Assemble the tridiagonal backward-Euler system for one closed-boundary bulk
salinity diffusion step.
"""
function assemble_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz,
            _assemble_column_salinity_system!,
            thermodynamics.auxiliary,
            thermodynamics.fields,
            grid,
            Δt)

    return nothing
end

"""
    solve_column_salinity_system!(thermodynamics)

Solve the assembled scalar salinity system into the bulk-salinity field.
"""
function solve_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics)
    solve!(thermodynamics.fields.bulk_salinity,
           thermodynamics.solvers.salinity_solver,
           thermodynamics.auxiliary.energy_rhs)

    return nothing
end

@kernel function _compute_column_metric_change!(scratch, field, grid)
    i, j, k = @index(Global, NTuple)
    metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
    moving_flux_difference =
        moving_face_displacement_flux(i, j, k+1, grid, field) -
        moving_face_displacement_flux(i, j, k, grid, field)
    moving_tendency = moving_flux_difference / Δzᶜᶜᶜ(i, j, k, grid)
    @inbounds scratch[i, j, k] = metric_ratio * field[i, j, k] + moving_tendency
end

@kernel function _copy_column_metric_change!(field, scratch)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = scratch[i, j, k]
end

function apply_column_metric_change!(field, scratch)
    grid = field.grid
    arch = architecture(grid)

    launch!(arch, grid, :xyz, _compute_column_metric_change!, scratch, field, grid)
    launch!(arch, grid, :xyz, _copy_column_metric_change!, field, scratch)

    return nothing
end

"""
    column_salinity_time_step!(thermodynamics, dt)

Advance prognostic bulk salinity by one scalar diffusion step. On a
`MutableVerticalDiscretization`, the step evolves layer-integrated salinity
with the moving vertical metric. The step is a no-op for prescribed salinity.
"""
column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics, Δt) =
    column_salinity_time_step!(thermodynamics,
                               thermodynamics.salinity_closure,
                               thermodynamics.salinity_transport,
                               Δt)

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrescribedBulkSalinity,
                                    salinity_transport,
                                    Δt)
    return nothing
end

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrognosticBulkSalinity,
                                    ::NoSalinityTransport,
                                    Δt)
    apply_column_metric_change!(thermodynamics.fields.bulk_salinity,
                                thermodynamics.auxiliary.energy_rhs)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    return nothing
end

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrognosticBulkSalinity,
                                    ::BrineSalinityDiffusion,
                                    Δt)
    throw(ArgumentError("BrineSalinityDiffusion is a marker for future brine-salinity transport; use BulkSalinityDiffusion for the implemented scalar bulk-salinity step."))
end

function column_salinity_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                    ::PrognosticBulkSalinity,
                                    ::BulkSalinityDiffusion,
                                    Δt)
    compute_column_salinity_diffusivity!(thermodynamics)
    assemble_column_salinity_system!(thermodynamics, Δt)
    solve_column_salinity_system!(thermodynamics)
    compute_column_thermodynamic_diagnostics!(thermodynamics)

    return nothing
end

@inline column_stefan_parameter(value::Number, i, j, k) = value
@inline column_stefan_parameter(value::AbstractField, i, j, k) = @inbounds value[i, j, k]

"""
    column_stefan_thickness_change(phase_transitions, sea_ice_density, residual_energy_flux, dt)

Return the conservative Stefan thickness change implied by a residual interface
energy flux. `residual_energy_flux` is positive for ice growth.
"""
@inline function column_stefan_thickness_change(phase_transitions,
                                                sea_ice_density,
                                                residual_energy_flux,
                                                Δt)
    volumetric_latent_heat = sea_ice_density * phase_transitions.reference_latent_heat
    return Δt * residual_energy_flux / volumetric_latent_heat
end

@kernel function _column_stefan_thickness_update!(ice_thickness,
                                                  phase_transitions,
                                                  sea_ice_density,
                                                  residual_energy_flux,
                                                  Δt)
    i, j = @index(Global, NTuple)

    @inbounds begin
        ρi = column_stefan_parameter(sea_ice_density, i, j, 1)
        δQ = column_stefan_parameter(residual_energy_flux, i, j, 1)
        ice_thickness[i, j, 1] += column_stefan_thickness_change(phase_transitions, ρi, δQ, Δt)
    end
end

"""
    column_stefan_thickness_update!(ice_thickness, phase_transitions, sea_ice_density, residual_energy_flux, dt)

Update `ice_thickness` conservatively from a residual interface energy flux.
The residual flux is positive for growth and negative for melt.
"""
function column_stefan_thickness_update!(ice_thickness,
                                         phase_transitions,
                                         sea_ice_density,
                                         residual_energy_flux,
                                         Δt)
    grid = ice_thickness.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy,
            _column_stefan_thickness_update!,
            ice_thickness,
            phase_transitions,
            sea_ice_density,
            residual_energy_flux,
            Δt)

    return nothing
end

"""
    column_stefan_thickness_budget(initial_thickness, final_thickness, phase_transitions, sea_ice_density, residual_energy_flux, dt)

Return a named tuple summarizing the conservative Stefan thickness update for
a scalar residual flux.
"""
function column_stefan_thickness_budget(initial_thickness,
                                        final_thickness,
                                        phase_transitions,
                                        sea_ice_density,
                                        residual_energy_flux,
                                        Δt)
    expected_change = column_stefan_thickness_change(phase_transitions,
                                                     sea_ice_density,
                                                     residual_energy_flux,
                                                     Δt)
    actual_change = final_thickness - initial_thickness
    residual = actual_change - expected_change

    return (initial = initial_thickness,
            final = final_thickness,
            actual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end

function check_column_remap_inputs(source_values, source_faces, target_faces)
    length(source_values) >= 1 ||
        throw(ArgumentError("source_values must contain at least one entry"))

    length(source_faces) == length(source_values) + 1 ||
        throw(ArgumentError("source_faces must have one more entry than source_values"))

    length(target_faces) >= 2 ||
        throw(ArgumentError("target_faces must contain at least two entries"))

    for k in 1:length(source_faces)-1
        source_faces[k+1] > source_faces[k] ||
            throw(ArgumentError("source_faces must be strictly increasing"))
    end

    for k in 1:length(target_faces)-1
        target_faces[k+1] > target_faces[k] ||
            throw(ArgumentError("target_faces must be strictly increasing"))
    end

    return nothing
end

"""
    conservative_column_remap!(target_values, source_values, source_faces, target_faces; fill_value=0)

Conservatively remap piecewise-constant layer averages from `source_faces` to
`target_faces`. Uncovered parts of the target column are filled with
`fill_value`, which represents newly exposed or newly grown ice.
"""
function conservative_column_remap!(target_values,
                                    source_values,
                                    source_faces,
                                    target_faces;
                                    fill_value = zero(eltype(target_values)))
    check_column_remap_inputs(source_values, source_faces, target_faces)

    length(target_values) == length(target_faces) - 1 ||
        throw(ArgumentError("target_values must have one fewer entry than target_faces"))

    for kt in eachindex(target_values)
        target_bottom = target_faces[kt]
        target_top = target_faces[kt+1]
        target_width = target_top - target_bottom
        integral = fill_value * target_width

        for ks in eachindex(source_values)
            source_bottom = source_faces[ks]
            source_top = source_faces[ks+1]
            overlap = min(target_top, source_top) - max(target_bottom, source_bottom)

            if overlap > 0
                integral += (source_values[ks] - fill_value) * overlap
            end
        end

        target_values[kt] = integral / target_width
    end

    return target_values
end

"""
    conservative_column_remap(source_values, source_faces, target_faces; fill_value=0)

Return layer averages on `target_faces` by conservative piecewise-constant
remapping from `source_faces`.
"""
function conservative_column_remap(source_values,
                                   source_faces,
                                   target_faces;
                                   fill_value = zero(eltype(source_values)))
    target_values = similar(source_values, length(target_faces) - 1)
    return conservative_column_remap!(target_values,
                                      source_values,
                                      source_faces,
                                      target_faces;
                                      fill_value)
end

@inline column_remap_parameter(value::Number, i, j) = value
@inline column_remap_parameter(value::AbstractField, i, j) = @inbounds value[i, j, 1]

"""
    column_energy_thickness_remap!(thermodynamics, source_faces, target_faces; fill_energy=0, bulk_salinity=nothing)

Conservatively remap layer-averaged internal energy in `thermodynamics` from
`source_faces` to `target_faces`, then recompute thermodynamic diagnostics.
Uncovered target intervals are filled with `fill_energy`, which can be a scalar
or a two-dimensional field. If `bulk_salinity` is provided, the bulk-salinity
field is reset with that value or profile before diagnostics are recomputed;
this is the fixed-salinity path used after thickness changes.
"""
function column_energy_thickness_remap!(thermodynamics::ColumnEnergyThermodynamics,
                                        source_faces,
                                        target_faces;
                                        fill_energy = zero(eltype(thermodynamics.fields.internal_energy.grid)),
                                        bulk_salinity = nothing)
    grid = thermodynamics.fields.internal_energy.grid
    on_cpu(grid) || error("column_energy_thickness_remap! currently supports CPU grids only.")
    check_column_remap_inputs(zeros(eltype(grid), size(grid, 3)), source_faces, target_faces)

    Nz = size(grid, 3)
    length(target_faces) == Nz + 1 ||
        throw(ArgumentError("target_faces must have one more entry than the number of vertical cells"))

    fields = thermodynamics.fields
    source_values = zeros(eltype(grid), Nz)
    target_values = zeros(eltype(grid), Nz)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        for k in 1:Nz
            @inbounds source_values[k] = fields.internal_energy[i, j, k]
        end

        fill_value = column_remap_parameter(fill_energy, i, j)
        conservative_column_remap!(target_values,
                                   source_values,
                                   source_faces,
                                   target_faces;
                                   fill_value)

        for k in 1:Nz
            @inbounds fields.internal_energy[i, j, k] = target_values[k]
        end
    end

    isnothing(bulk_salinity) || set_column_field!(fields.bulk_salinity, bulk_salinity)
    compute_column_thermodynamic_diagnostics!(thermodynamics)

    return nothing
end

"""
    column_layer_integral(values, faces)

Return the integral of piecewise-constant layer averages over `faces`.
"""
function column_layer_integral(values, faces)
    length(values) >= 1 ||
        throw(ArgumentError("values must contain at least one entry"))

    length(faces) == length(values) + 1 ||
        throw(ArgumentError("faces must have one more entry than values"))

    FT = promote_type(eltype(values), eltype(faces))
    integral = zero(FT)

    for k in eachindex(values)
        width = faces[k+1] - faces[k]
        width > 0 ||
            throw(ArgumentError("faces must be strictly increasing"))
        integral += values[k] * width
    end

    return integral
end

@inline column_basal_stefan_flux(boundary, fields, auxiliary, grid, relation, i, j) = zero(eltype(grid))

@inline function column_basal_stefan_flux(boundary::ColumnDirichletBoundary, fields, auxiliary, grid, relation, i, j)
    conductance = column_boundary_temperature_conductance(i, j, 1, 1, grid, auxiliary)
    T₁ = @inbounds fields.temperature[i, j, 1]
    Tᵇ = column_dirichlet_temperature(boundary, i, j, grid, relation)
    return conductance * (Tᵇ - T₁)
end

# Unconsolidated top heat loss (upward-positive), mirroring SlabThermodynamics so a thin column behaves like a
# slab. A direct-flux/melting top exchanges its external flux; a Dirichlet (prescribed-temperature or ocean) top
# additionally loses heat by single-layer conduction k/h (Tᵇ − Tᵗ) — the same conduction the slab wraps into its
# top flux — so a thin conduction-driven column still grows instead of stalling.
@inline column_unconsolidated_top_flux(bc, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields) =
    getflux(ext, i, j, grid, Tᵇ, clock, model_fields)

@inline function column_unconsolidated_top_flux(bc::ColumnDirichletBoundary, ext, conductivity, Tᵇ, h, i, j, grid, fields, relation, clock, model_fields)
    Tᵗ = column_dirichlet_temperature(bc, i, j, grid, relation)
    S = @inbounds fields.bulk_salinity[i, j, 1]
    k = ice_thermal_conductivity(conductivity, Tᵗ, S)
    return k / h * (Tᵇ - Tᵗ) + getflux(ext, i, j, grid, Tᵗ, clock, model_fields)
end

@kernel function _column_stefan_volume_update!(ice_thickness, ice_concentration, consolidation_thickness,
                                               sea_ice_density, phase_transitions, z, fields, auxiliary, grid,
                                               heat_boundary_conditions, external_heat_fluxes, relation, conductivity,
                                               clock, model_fields, Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        hⁿ = ice_thickness[i, j, 1]
        ℵⁿ = ice_concentration[i, j, 1]
        hᶜ = consolidation_thickness[i, j, 1]
        ρᵢ = column_stefan_parameter(sea_ice_density, i, j, 1)

        # Roll the interface heights into the previous slot for the next moving-grid solve.
        z.hs⁻[i, j, 1] = z.hsⁿ[i, j, 1]
        z.hb⁻[i, j, 1] = z.hbⁿ[i, j, 1]

        if hⁿ ≥ hᶜ
            # Consolidated: resolved conduction drives basal congelation and the surface melt residual.
            basal_flux = column_basal_stefan_flux(heat_boundary_conditions.bottom, fields, auxiliary, grid, relation, i, j)
            surface_flux = column_surface_stefan_residual_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                                               i, j, Nz, grid, auxiliary, fields, relation, clock, model_fields, Δt)

            δhᵇ = column_stefan_thickness_change(phase_transitions, ρᵢ, basal_flux, Δt)
            δhˢ = column_stefan_thickness_change(phase_transitions, ρᵢ, surface_flux, Δt)

            # Per-ice thickness velocity (no ℵ premultiplication), matching SlabThermodynamics — `ice_volume_update`
            # and the concentration step grow ℵ from open water, so a column nucleates ice from ℵ = 0 like the slab.
            ∂t_V = (δhᵇ + δhˢ) / Δt
            surface_shift = δhˢ
        else
            # Unconsolidated: no resolved interior profile, so the column behaves as a single-layer slab. Net growth
            # is ℰ⁻¹ (Qᵘ − Qᵇ) with the shared upward-positive fluxes, where a Dirichlet top contributes its
            # single-layer conduction k/h just as SlabThermodynamics does.
            Tᵇ = column_dirichlet_temperature(heat_boundary_conditions.bottom, i, j, grid, relation)
            ℰ = ρᵢ * latent_heat(phase_transitions, Tᵇ)
            Qᵘ = column_unconsolidated_top_flux(heat_boundary_conditions.top, external_heat_fluxes.top,
                                                conductivity, Tᵇ, hⁿ, i, j, grid, fields, relation, clock, model_fields)
            Qᵇ = getflux(external_heat_fluxes.bottom, i, j, grid, Tᵇ, clock, model_fields)

            ∂t_V = (Qᵘ - Qᵇ) / ℰ
            surface_shift = zero(eltype(grid))
        end

        h⁺, ℵ⁺ = ice_volume_update(ProportionalEvolution(), ∂t_V, hⁿ, ℵⁿ, hᶜ, Δt)

        # Keep the grid metric consistent with the (area-averaged) ice thickness: surface melt moves the surface,
        # and the base is placed so the column height equals the updated thickness. For ℵ = 1 this is exactly the
        # basal Stefan move (z.hbⁿ -= δhᵇ); for evolving ℵ it follows the redistributed thickness from ice_volume_update.
        z.hsⁿ[i, j, 1] = z.hs⁻[i, j, 1] + surface_shift
        z.hbⁿ[i, j, 1] = z.hsⁿ[i, j, 1] - h⁺

        ice_thickness[i, j, 1] = h⁺
        ice_concentration[i, j, 1] = ℵ⁺
    end
end

function column_stefan_volume_update!(thermodynamics::ColumnEnergyThermodynamics, model, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    launch!(architecture(grid), grid, :xy, _column_stefan_volume_update!,
            model.ice_thickness, model.ice_concentration, model.ice_consolidation_thickness,
            model.sea_ice_density, model.phase_transitions, sea_ice_discretization(grid),
            thermodynamics.fields, thermodynamics.auxiliary, grid,
            thermodynamics.heat_boundary_conditions, model.external_heat_fluxes, thermodynamics.relation,
            thermal_conductivity(thermodynamics.energy_transport),
            model.clock, fields(model), Δt)
    return nothing
end

# The surface energy solve is a single tridiagonal pass for direct-flux and Dirichlet tops; a
# `MeltingConstrainedFluxBalance` top adds the outer surface-temperature iteration around it.
function column_surface_energy_solve!(thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    column_surface_energy_solve!(thermodynamics.heat_boundary_conditions.top,
                                 thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    return nothing
end

function column_surface_energy_solve!(top_boundary, thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    solve_column_energy_system!(thermodynamics)
    return nothing
end

function column_surface_energy_solve!(top_boundary::MeltingConstrainedFluxBalance,
                                      thermodynamics::ColumnEnergyThermodynamics,
                                      external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt;
                                      max_iterations = 20)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)
    tolerance = sqrt(eps(eltype(grid)))
    surface_temperature = thermodynamics.auxiliary.surface_temperature
    internal_energy = thermodynamics.fields.internal_energy

    # Each outer pass re-solves from the same start-of-step enthalpy with the latest surface coupling, exactly as
    # Icepack `temperature_changes` rebuilds the system from the initial temperature every iteration. Restoring it
    # (rather than accumulating onto the previous solve) is what keeps the lagged top flux applied once per step.
    start_of_step_energy = copy(interior(internal_energy))

    for iteration in 1:max_iterations
        previous = Array(interior(surface_temperature))
        interior(internal_energy) .= start_of_step_energy
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        assemble_column_energy_system!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
        solve_column_energy_system!(thermodynamics)
        compute_column_thermodynamic_diagnostics!(thermodynamics)
        launch!(arch, grid, :xy, _update_column_surface_temperature!,
                thermodynamics.auxiliary, thermodynamics.fields, grid,
                external_heat_fluxes.top, thermodynamics.relation, clock, model_fields)
        maximum(abs.(Array(interior(surface_temperature)) .- previous)) <= tolerance && break
    end

    return nothing
end

# Standalone solves treat every column as consolidated (`nothing`); the coupled step passes the model's
# consolidation thickness so thin columns take the slab-balance regime.
column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics, external_heat_fluxes, clock, model_fields, Δt) =
    column_energy_time_step!(thermodynamics, external_heat_fluxes, clock, model_fields, nothing, Δt)

function column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics,
                                  external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    column_surface_energy_solve!(thermodynamics, external_heat_fluxes, clock, model_fields, consolidation_thickness, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    column_salinity_time_step!(thermodynamics, Δt)

    return nothing
end

function thermodynamic_time_step!(model,
                                  thermodynamics::ColumnEnergyThermodynamics,
                                  ::Nothing,
                                  Δt)
    column_energy_time_step!(thermodynamics, model.external_heat_fluxes, model.clock, fields(model),
                             model.ice_consolidation_thickness, Δt)
    column_stefan_volume_update!(thermodynamics, model, Δt)
    return nothing
end

"""
    column_integrated_energy(thermodynamics)

Return the vertical integral of volumetric internal energy over every column
cell in `thermodynamics`.
"""
function column_integrated_energy(thermodynamics::ColumnEnergyThermodynamics)
    E = thermodynamics.fields.internal_energy
    grid = E.grid
    total_energy = zero(eltype(grid))

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        total_energy += @inbounds E[i, j, k] * Δzᶜᶜᶜ(i, j, k, grid)
    end

    return total_energy
end

"""
    column_integrated_salinity(thermodynamics)

Return the vertical integral of bulk salinity over every column cell in
`thermodynamics`.
"""
function column_integrated_salinity(thermodynamics::ColumnEnergyThermodynamics)
    S = thermodynamics.fields.bulk_salinity
    grid = S.grid
    total_salinity = zero(eltype(grid))

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        total_salinity += @inbounds S[i, j, k] * Δzᶜᶜᶜ(i, j, k, grid)
    end

    return total_salinity
end

function column_boundary_energy_flux_difference(thermodynamics::ColumnEnergyThermodynamics,
                                                external_heat_fluxes, clock, model_fields)
    grid = thermodynamics.fields.internal_energy.grid
    temperature = thermodynamics.fields.temperature
    Nz = size(grid, 3)
    difference = zero(eltype(grid))

    # Upward-positive fluxes: net energy gained by the column is (ocean-in) − (top-out) = Qᵇ − Qᵘ.
    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        T_top = @inbounds temperature[i, j, Nz]
        T_bottom = @inbounds temperature[i, j, 1]
        top = getflux(external_heat_fluxes.top, i, j, grid, T_top, clock, model_fields)
        bottom = getflux(external_heat_fluxes.bottom, i, j, grid, T_bottom, clock, model_fields)
        difference += bottom - top
    end

    return difference
end

function column_shortwave_flux_difference(thermodynamics::ColumnEnergyThermodynamics)
    compute_column_shortwave_flux!(thermodynamics)

    I = thermodynamics.auxiliary.shortwave_flux
    grid = thermodynamics.fields.internal_energy.grid
    Nz = size(grid, 3)
    difference = zero(eltype(grid))

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        difference += @inbounds I[i, j, Nz+1] - I[i, j, 1]
    end

    return difference
end

@inline function budget_relative_residual(residual, actual_change, expected_change)
    scale = max(abs(actual_change), abs(expected_change), one(abs(residual)))
    return abs(residual) / scale
end

"""
    column_energy_budget(thermodynamics, external_heat_fluxes, clock, model_fields, initial_energy, dt;
                         final_energy = column_integrated_energy(thermodynamics),
                         surface_stefan_residual_flux = 0)

Return a named tuple summarizing the column-integrated energy budget over one
step with duration `dt`. The boundary flux change is the requested (uncapped)
`external_heat_fluxes` difference; surface melt at a `MeltingConstrainedFluxBalance`
boundary is accounted separately through `surface_stefan_residual_flux`.
"""
function column_energy_budget(thermodynamics::ColumnEnergyThermodynamics,
                              external_heat_fluxes, clock, model_fields,
                              initial_energy,
                              Δt;
                              final_energy = column_integrated_energy(thermodynamics),
                              surface_stefan_residual_flux = 0)
    boundary_flux_change = Δt * column_boundary_energy_flux_difference(thermodynamics, external_heat_fluxes, clock, model_fields)
    shortwave_flux_change = Δt * column_shortwave_flux_difference(thermodynamics)
    surface_stefan_residual_change = Δt * surface_stefan_residual_flux
    expected_change = boundary_flux_change + shortwave_flux_change + surface_stefan_residual_change
    actual_change = final_energy - initial_energy
    residual = actual_change - expected_change

    return (initial = initial_energy,
            final = final_energy,
            actual_change,
            boundary_flux_change,
            shortwave_flux_change,
            surface_stefan_residual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end

"""
    column_salt_budget(thermodynamics, initial_salt, dt; final_salt = column_integrated_salinity(thermodynamics))

Return a named tuple summarizing the closed-boundary column-integrated salt
budget over one step with duration `dt`.
"""
function column_salt_budget(thermodynamics::ColumnEnergyThermodynamics,
                            initial_salt,
                            Δt;
                            final_salt = column_integrated_salinity(thermodynamics))
    expected_change = zero(final_salt - initial_salt)
    actual_change = final_salt - initial_salt
    residual = actual_change - expected_change

    return (initial = initial_salt,
            final = final_salt,
            actual_change,
            expected_change,
            residual,
            relative_residual = budget_relative_residual(residual,
                                                         actual_change,
                                                         expected_change))
end

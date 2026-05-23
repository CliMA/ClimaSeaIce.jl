import Oceananigans: prognostic_state, restore_prognostic_state!
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: CPU, architecture
using Oceananigans.Fields: AbstractField
using Oceananigans.Grids: ZDirection, rnode, znode
using Oceananigans.Operators: Δzᵃᵃᶜ, Δzᵃᵃᶠ, σ⁻, σⁿ
using Oceananigans.Solvers: BatchedTridiagonalSolver, solve!
using Oceananigans.Utils: launch!

@inline on_cpu(grid) = architecture(grid) isa CPU

@inline function current_column_cell_thickness(i, j, k, grid)
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Center())
    return Δzᵃᵃᶜ(i, j, k, grid) * σᶜᶜⁿ
end

@inline function current_column_face_spacing(i, j, k, grid)
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Face())
    return Δzᵃᵃᶠ(i, j, k, grid) * σᶜᶜⁿ
end

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

"""
    InsulatingBoundary()

Column energy boundary condition with zero normal energy flux.
"""
struct InsulatingBoundary end

"""
    PrescribedEnergyFlux([FT=Oceananigans.defaults.FloatType; flux=0])

Column energy boundary condition with prescribed face energy flux. Flux is
positive in the increasing vertical-coordinate direction.
"""
struct PrescribedEnergyFlux{F}
    flux :: F
end

PrescribedEnergyFlux(FT::DataType=Oceananigans.defaults.FloatType; flux = 0) =
    PrescribedEnergyFlux(convert(FT, flux))

"""
    PrescribedEnergyFluxBoundaryEnergy([FT=Oceananigans.defaults.FloatType; flux=0, boundary_energy=0])

Column energy boundary condition with a prescribed face energy flux and a
prescribed volumetric internal energy for material swept into the column by a
moving boundary. This reduces to [`PrescribedEnergyFlux`](@ref) on static grids.
"""
struct PrescribedEnergyFluxBoundaryEnergy{F, E}
    flux :: F
    boundary_energy :: E
end

PrescribedEnergyFluxBoundaryEnergy(FT::DataType=Oceananigans.defaults.FloatType;
                                   flux = 0,
                                   boundary_energy = 0) =
    PrescribedEnergyFluxBoundaryEnergy(convert(FT, flux),
                                       convert(FT, boundary_energy))

function Adapt.adapt_structure(to, boundary::PrescribedEnergyFluxBoundaryEnergy)
    return PrescribedEnergyFluxBoundaryEnergy(Adapt.adapt(to, boundary.flux),
                                             Adapt.adapt(to, boundary.boundary_energy))
end

"""
    MeltingLimitedSurfaceFlux([FT=Oceananigans.defaults.FloatType; flux=0])

Top energy boundary condition with an imposed surface flux that is capped by
the complete-melt energy of the top cell. The applied flux warms the column
column up to the complete-melt threshold; the excess is available as a
conservative Stefan residual for surface melt.
"""
struct MeltingLimitedSurfaceFlux{F}
    flux :: F
end

MeltingLimitedSurfaceFlux(FT::DataType=Oceananigans.defaults.FloatType; flux = 0) =
    MeltingLimitedSurfaceFlux(convert(FT, flux))

"""
    ColumnBoundaryConditions(; top=InsulatingBoundary(), bottom=InsulatingBoundary())

Top and bottom energy boundary conditions for the column energy
solve.
"""
struct ColumnBoundaryConditions{T, B}
    top :: T
    bottom :: B
end

ColumnBoundaryConditions(; top = InsulatingBoundary(),
                         bottom = InsulatingBoundary()) =
    ColumnBoundaryConditions(top, bottom)

struct ColumnThermodynamicFields{E, S, T, LF, BS, TE, TS}
    internal_energy :: E
    bulk_salinity :: S
    temperature :: T
    liquid_fraction :: LF
    brine_salinity :: BS
    temperature_energy_derivative :: TE
    temperature_salinity_derivative :: TS
end

struct ColumnAuxiliaryFields{K, KE, KS, DE, CS, I, RHS, A, B, C}
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
    boundary_conditions :: BC
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
                                 Adapt.adapt(to, auxiliary.upper_diagonal))
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
                                      Adapt.adapt(to, thermodynamics.boundary_conditions),
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

    return ColumnAuxiliaryFields(thermal_conductivity,
                                 energy_diffusivity,
                                 salinity_diffusivity,
                                 effective_energy_diffusivity,
                                 salinity_coupling_diffusivity,
                                 shortwave_flux,
                                 energy_rhs,
                                 lower_diagonal,
                                 diagonal,
                                 upper_diagonal)
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

function ColumnEnergyThermodynamics(grid;
                                    relation = QuadraticLiquidusEnergyRelation(eltype(grid)),
                                    salinity_closure = PrognosticBulkSalinity(),
                                    energy_transport = ConductiveTemperatureTransport(eltype(grid)),
                                    salinity_transport = NoSalinityTransport(),
                                    shortwave_absorption = NoShortwaveAbsorption(),
                                    boundary_conditions = ColumnBoundaryConditions(),
                                    fields = column_thermodynamic_fields(grid),
                                    auxiliary = column_auxiliary_fields(grid),
                                    solvers = column_solvers(grid, auxiliary))
    initialize_salinity!(fields, salinity_closure)

    return ColumnEnergyThermodynamics(relation,
                                      salinity_closure,
                                      energy_transport,
                                      salinity_transport,
                                      shortwave_absorption,
                                      boundary_conditions,
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
                                                     boundary_conditions = ColumnBoundaryConditions(),
                                                     kw...)
    return ColumnEnergyThermodynamics(grid;
                                      relation,
                                      salinity_closure = PrescribedBulkSalinity(salinity_profile),
                                      energy_transport,
                                      salinity_transport = NoSalinityTransport(),
                                      shortwave_absorption,
                                      boundary_conditions,
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
                                                boundary_conditions = ColumnBoundaryConditions(),
                                                kw...)
    return ColumnEnergyThermodynamics(grid;
                                      relation,
                                      salinity_closure = PrognosticBulkSalinity(),
                                      energy_transport,
                                      salinity_transport,
                                      shortwave_absorption,
                                      boundary_conditions,
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

function _compute_column_internal_energy_cpu!(fields, relation)
    grid = fields.internal_energy.grid

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds fields.internal_energy[i, j, k] =
            internal_energy(relation,
                            fields.temperature[i, j, k],
                            fields.bulk_salinity[i, j, k])
    end

    return nothing
end

"""
    compute_column_internal_energy!(thermodynamics)

Fill the internal-energy field from the current temperature and bulk-salinity
fields.
"""
function compute_column_internal_energy!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_internal_energy_cpu!(thermodynamics.fields,
                                             thermodynamics.relation)
    else
        launch!(arch, grid, :xyz,
                _compute_column_internal_energy!,
                thermodynamics.fields,
                thermodynamics.relation)
    end

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

function _compute_column_thermodynamic_diagnostics_cpu!(fields, relation)
    grid = fields.internal_energy.grid

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
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

    return nothing
end

"""
    compute_column_thermodynamic_diagnostics!(thermodynamics)

Fill temperature, liquid fraction, brine salinity, and thermodynamic derivative
diagnostics from internal energy and bulk salinity.
"""
function compute_column_thermodynamic_diagnostics!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_thermodynamic_diagnostics_cpu!(thermodynamics.fields,
                                                       thermodynamics.relation)
    else
        launch!(arch, grid, :xyz,
                _compute_column_thermodynamic_diagnostics!,
                thermodynamics.fields,
                thermodynamics.relation)
    end

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
                                          current_column_cell_thickness(i, j, k-1, grid),
                                          T[i, j, k],
                                          S[i, j, k],
                                          current_column_cell_thickness(i, j, k, grid))

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

function _compute_column_transport_coefficients_cpu!(auxiliary, fields, energy_transport)
    grid = fields.internal_energy.grid
    Nz = size(grid, 3)

    kappa_T = thermal_conductivity(energy_transport)
    kappa_E = energy_diffusivity(energy_transport)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            T = fields.temperature
            S = fields.bulk_salinity

            boundary_conductivity = ice_thermal_conductivity(kappa_T, T[i, j, 1], S[i, j, 1])
            TE = fields.temperature_energy_derivative[i, j, 1]
            TS = fields.temperature_salinity_derivative[i, j, 1]

            auxiliary.thermal_conductivity[i, j, 1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, 1] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, 1] = boundary_conductivity * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, 1] = boundary_conductivity * TS

            for k in 2:Nz
                TE = (fields.temperature_energy_derivative[i, j, k-1] +
                      fields.temperature_energy_derivative[i, j, k]) / 2

                TS = (fields.temperature_salinity_derivative[i, j, k-1] +
                      fields.temperature_salinity_derivative[i, j, k]) / 2

                face_conductivity =
                    face_thermal_conductivity(kappa_T,
                                              T[i, j, k-1],
                                              S[i, j, k-1],
                                              current_column_cell_thickness(i, j, k-1, grid),
                                              T[i, j, k],
                                              S[i, j, k],
                                              current_column_cell_thickness(i, j, k, grid))

                auxiliary.thermal_conductivity[i, j, k] = face_conductivity
                auxiliary.energy_diffusivity[i, j, k] = kappa_E
                auxiliary.effective_energy_diffusivity[i, j, k] = face_conductivity * TE + kappa_E
                auxiliary.salinity_coupling_diffusivity[i, j, k] = face_conductivity * TS
            end

            boundary_conductivity = ice_thermal_conductivity(kappa_T, T[i, j, Nz], S[i, j, Nz])
            TE = fields.temperature_energy_derivative[i, j, Nz]
            TS = fields.temperature_salinity_derivative[i, j, Nz]

            auxiliary.thermal_conductivity[i, j, Nz+1] = boundary_conductivity
            auxiliary.energy_diffusivity[i, j, Nz+1] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, Nz+1] = boundary_conductivity * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, Nz+1] = boundary_conductivity * TS
        end
    end

    return nothing
end

"""
    compute_column_transport_coefficients!(thermodynamics)

Compute face-centered energy transport coefficients for the semi-implicit
energy solve.
"""
function compute_column_transport_coefficients!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_transport_coefficients_cpu!(thermodynamics.auxiliary,
                                                    thermodynamics.fields,
                                                    thermodynamics.energy_transport)
    else
        launch!(arch, grid, :xyz,
                _compute_column_transport_coefficients!,
                thermodynamics.auxiliary,
                thermodynamics.fields,
                grid,
                thermodynamics.energy_transport)
    end

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

function _compute_column_salinity_diffusivity_cpu!(auxiliary, fields, salinity_transport)
    grid = fields.bulk_salinity.grid
    Nz = size(grid, 3)
    kappa_S = salinity_diffusivity(salinity_transport)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            auxiliary.salinity_diffusivity[i, j, 1] = zero(kappa_S)

            for k in 2:Nz
                auxiliary.salinity_diffusivity[i, j, k] = kappa_S
            end

            auxiliary.salinity_diffusivity[i, j, Nz+1] = zero(kappa_S)
        end
    end

    return nothing
end

"""
    compute_column_salinity_diffusivity!(thermodynamics)

Compute face-centered bulk-salinity diffusivity for the scalar salinity solve.
"""
function compute_column_salinity_diffusivity!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_salinity_diffusivity_cpu!(thermodynamics.auxiliary,
                                                  thermodynamics.fields,
                                                  thermodynamics.salinity_transport)
    else
        launch!(arch, grid, :xyz,
                _compute_column_salinity_diffusivity!,
                thermodynamics.auxiliary,
                thermodynamics.fields,
                thermodynamics.salinity_transport)
    end

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

function _compute_column_shortwave_flux_cpu!(auxiliary, grid, shortwave_absorption)
    Nz = size(grid, 3)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        z_top = znode(i, j, Nz+1, grid, Center(), Center(), Face())

        @inbounds begin
            z_bottom = znode(i, j, 1, grid, Center(), Center(), Face())
            auxiliary.shortwave_flux[i, j, 1] =
                shortwave_flux(shortwave_absorption, z_bottom, z_top)

            for k in 2:Nz
                z_face = znode(i, j, k, grid, Center(), Center(), Face())
                auxiliary.shortwave_flux[i, j, k] =
                    shortwave_flux(shortwave_absorption, z_face, z_top)
            end

            auxiliary.shortwave_flux[i, j, Nz+1] =
                shortwave_flux(shortwave_absorption, z_top, z_top)
        end
    end

    return nothing
end

"""
    compute_column_shortwave_flux!(thermodynamics)

Compute face-centered shortwave fluxes for the current shortwave absorption
closure.
"""
function compute_column_shortwave_flux!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_shortwave_flux_cpu!(thermodynamics.auxiliary,
                                            grid,
                                            thermodynamics.shortwave_absorption)
    else
        launch!(arch, grid, :xyz,
                _compute_column_shortwave_flux!,
                thermodynamics.auxiliary,
                grid,
                thermodynamics.shortwave_absorption)
    end

    return nothing
end

@kernel function _compute_column_surface_stefan_residual_flux!(residual_flux,
                                                               fields,
                                                               grid,
                                                               boundary_conditions,
                                                               relation,
                                                               Δt)
    i, j = @index(Global, NTuple)
    Nz = size(grid, 3)

    @inbounds begin
        E = fields.internal_energy[i, j, Nz]
        S = fields.bulk_salinity[i, j, Nz]
        Δz = current_column_cell_thickness(i, j, Nz, grid)

        residual_flux[i, j, 1] =
            column_surface_stefan_residual_flux(boundary_conditions.top,
                                                relation,
                                                E,
                                                S,
                                                Δz,
                                                Δt)
    end
end

function _compute_column_surface_stefan_residual_flux_cpu!(residual_flux,
                                                           fields,
                                                           grid,
                                                           boundary_conditions,
                                                           relation,
                                                           Δt)
    Nz = size(grid, 3)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            E = fields.internal_energy[i, j, Nz]
            S = fields.bulk_salinity[i, j, Nz]
            Δz = current_column_cell_thickness(i, j, Nz, grid)

            residual_flux[i, j, 1] =
                column_surface_stefan_residual_flux(boundary_conditions.top,
                                                    relation,
                                                    E,
                                                    S,
                                                    Δz,
                                                    Δt)
        end
    end

    return nothing
end

"""
    compute_column_surface_stefan_residual_flux!(residual_flux, thermodynamics, dt)

Fill `residual_flux` with the top-surface Stefan residual implied by a
melting-limited surface boundary. The residual is positive for growth and
negative for melt, matching [`column_stefan_thickness_update!`](@ref).
"""
function compute_column_surface_stefan_residual_flux!(residual_flux,
                                                      thermodynamics::ColumnEnergyThermodynamics,
                                                      Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _compute_column_surface_stefan_residual_flux_cpu!(residual_flux,
                                                          thermodynamics.fields,
                                                          grid,
                                                          thermodynamics.boundary_conditions,
                                                          thermodynamics.relation,
                                                          Δt)
    else
        launch!(arch, grid, :xy,
                _compute_column_surface_stefan_residual_flux!,
                residual_flux,
                thermodynamics.fields,
                grid,
                thermodynamics.boundary_conditions,
                thermodynamics.relation,
                Δt)
    end

    return nothing
end

"""
    column_surface_stefan_residual_flux(thermodynamics, dt)

Return the column-summed top-surface Stefan residual implied by the current
state and top boundary condition.
"""
function column_surface_stefan_residual_flux(thermodynamics::ColumnEnergyThermodynamics, Δt)
    fields = thermodynamics.fields
    grid = fields.internal_energy.grid
    top = thermodynamics.boundary_conditions.top
    residual = zero(eltype(grid))
    Nz = size(grid, 3)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            E = fields.internal_energy[i, j, Nz]
            S = fields.bulk_salinity[i, j, Nz]
            Δz = current_column_cell_thickness(i, j, Nz, grid)

            residual += column_surface_stefan_residual_flux(top,
                                                            thermodynamics.relation,
                                                            E,
                                                            S,
                                                            Δz,
                                                            Δt)
        end
    end

    return residual
end

@inline function diffusion_factor(i, j, kc, kf, grid, diffusivity, Δt)
    return Δt * @inbounds(diffusivity[i, j, kf]) /
           (current_column_cell_thickness(i, j, kc, grid) *
            current_column_face_spacing(i, j, kf, grid))
end

@inline function previous_to_current_column_metric_ratio(i, j, k, grid)
    σᶜᶜ⁻ = σ⁻(i, j, k, grid, Center(), Center(), Center())
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Center())
    return σᶜᶜ⁻ / σᶜᶜⁿ
end

@inline function column_face_displacement(i, j, k, grid)
    σᶜᶜ⁻ = σ⁻(i, j, k, grid, Center(), Center(), Face())
    σᶜᶜⁿ = σⁿ(i, j, k, grid, Center(), Center(), Face())
    r = rnode(i, j, k, grid, Center(), Center(), Face())
    return (σᶜᶜⁿ - σᶜᶜ⁻) * r
end

@inline moving_boundary_value(value::Number, i, j) = value
@inline moving_boundary_value(value, i, j) = @inbounds value[i, j]

@inline function moving_bottom_boundary_value(boundary, i, j, default)
    return default
end

@inline function moving_top_boundary_value(boundary, i, j, default)
    return default
end

@inline function moving_bottom_boundary_value(boundary::PrescribedEnergyFluxBoundaryEnergy,
                                              i, j, default)
    return moving_boundary_value(boundary.boundary_energy, i, j)
end

@inline function moving_top_boundary_value(boundary::PrescribedEnergyFluxBoundaryEnergy,
                                           i, j, default)
    return moving_boundary_value(boundary.boundary_energy, i, j)
end

@inline function moving_face_value(i, j, k, grid, field, face_displacement,
                                   bottom_boundary_value,
                                   top_boundary_value)
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

@inline function moving_face_displacement_flux(i, j, k, grid, field,
                                               bottom_boundary,
                                               top_boundary)
    displacement = column_face_displacement(i, j, k, grid)
    Nz = size(grid, 3)
    bottom_default = @inbounds field[i, j, 1]
    top_default = @inbounds field[i, j, Nz]
    bottom_boundary_value = moving_bottom_boundary_value(bottom_boundary, i, j, bottom_default)
    top_boundary_value = moving_top_boundary_value(top_boundary, i, j, top_default)
    value = moving_face_value(i, j, k, grid, field, displacement,
                              bottom_boundary_value, top_boundary_value)
    return displacement * value
end

@inline function salinity_coupling_flux(i, j, k, grid, C, S)
    return @inbounds(C[i, j, k] * (S[i, j, k] - S[i, j, k-1]) /
                     current_column_face_spacing(i, j, k, grid))
end

@inline column_boundary_energy_flux(::InsulatingBoundary) = 0
@inline column_boundary_energy_flux(boundary::PrescribedEnergyFlux) = boundary.flux
@inline column_boundary_energy_flux(boundary::PrescribedEnergyFluxBoundaryEnergy) = boundary.flux
@inline column_boundary_energy_flux(boundary::MeltingLimitedSurfaceFlux) = boundary.flux

@inline column_requested_surface_energy_flux(boundary) = column_boundary_energy_flux(boundary)

@inline function column_boundary_temperature_conductance(i, j, kf, kc, grid, auxiliary)
    K = @inbounds auxiliary.thermal_conductivity[i, j, kf]
    half_cell_thickness = current_column_cell_thickness(i, j, kc, grid) / 2
    return K / half_cell_thickness
end

@inline column_bottom_boundary_energy_factor(boundary, i, j, k, grid, auxiliary, fields, Δt) = zero(Δt)
@inline column_top_boundary_energy_factor(boundary, i, j, k, grid, auxiliary, fields, Δt) = zero(Δt)

@inline function column_bottom_boundary_energy_factor(boundary::PrescribedTemperature,
                                                      i, j, k, grid, auxiliary, fields, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / current_column_cell_thickness(i, j, k, grid)
end

@inline function column_top_boundary_energy_factor(boundary::PrescribedTemperature,
                                                   i, j, k, grid, auxiliary, fields, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    return Δt * conductance * TE / current_column_cell_thickness(i, j, k, grid)
end

@inline function column_bottom_boundary_energy_flux(boundary,
                                                    i, j, k, grid, auxiliary, fields, relation, Δt)
    return column_boundary_energy_flux(boundary, i, j, k, grid, fields, relation, Δt)
end

@inline function column_top_boundary_energy_flux(boundary,
                                                 i, j, k, grid, auxiliary, fields, relation, Δt)
    return column_boundary_energy_flux(boundary, i, j, k, grid, fields, relation, Δt)
end

@inline function column_bottom_boundary_energy_flux(boundary::PrescribedTemperature,
                                                    i, j, k, grid, auxiliary, fields, relation, Δt)
    conductance = column_boundary_temperature_conductance(i, j, 1, k, grid, auxiliary)
    T = @inbounds fields.temperature[i, j, k]
    E = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    boundary_temperature = column_prescribed_temperature(boundary, i, j)
    return conductance * (T - boundary_temperature - TE * E)
end

@inline function column_top_boundary_energy_flux(boundary::PrescribedTemperature,
                                                 i, j, k, grid, auxiliary, fields, relation, Δt)
    Nz = size(grid, 3)
    conductance = column_boundary_temperature_conductance(i, j, Nz+1, k, grid, auxiliary)
    T = @inbounds fields.temperature[i, j, k]
    E = @inbounds fields.internal_energy[i, j, k]
    TE = @inbounds fields.temperature_energy_derivative[i, j, k]
    boundary_temperature = column_prescribed_temperature(boundary, i, j)
    return conductance * (boundary_temperature - T + TE * E)
end

@inline function column_energy_to_complete_melt(relation, E, S, Δz)
    return (complete_melt_energy(relation, S) - E) * Δz
end

@inline function column_surface_energy_flux(boundary::MeltingLimitedSurfaceFlux,
                                            relation,
                                            E,
                                            S,
                                            Δz,
                                            Δt)
    requested_flux = column_requested_surface_energy_flux(boundary)
    available_flux = column_energy_to_complete_melt(relation, E, S, Δz) / Δt
    nonnegative_available_flux = max(available_flux, zero(available_flux))
    return min(requested_flux, nonnegative_available_flux)
end

@inline column_surface_energy_flux(boundary, relation, E, S, Δz, Δt) =
    column_requested_surface_energy_flux(boundary)

@inline function column_surface_stefan_residual_flux(boundary::MeltingLimitedSurfaceFlux,
                                                     relation,
                                                     E,
                                                     S,
                                                     Δz,
                                                     Δt)
    return column_surface_energy_flux(boundary, relation, E, S, Δz, Δt) -
           column_requested_surface_energy_flux(boundary)
end

@inline column_surface_stefan_residual_flux(boundary, relation, E, S, Δz, Δt) = zero(E)

@inline function column_boundary_energy_flux(boundary, i, j, k, grid, fields, relation, Δt)
    Δz = current_column_cell_thickness(i, j, k, grid)
    E = @inbounds fields.internal_energy[i, j, k]
    S = @inbounds fields.bulk_salinity[i, j, k]
    return column_surface_energy_flux(boundary, relation, E, S, Δz, Δt)
end

@kernel function _assemble_column_energy_system!(auxiliary, fields, grid, boundary_conditions, relation, Δt)
    i, j, k = @index(Global, NTuple)
    Nz = size(grid, 3)

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
            column_bottom_boundary_energy_factor(boundary_conditions.bottom,
                                                 i, j, k, grid, auxiliary, fields, Δt)
        else
            zero(Δt)
        end

        top_factor = if k == Nz
            column_top_boundary_energy_factor(boundary_conditions.top,
                                              i, j, k, grid, auxiliary, fields, Δt)
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
            column_bottom_boundary_energy_flux(boundary_conditions.bottom,
                                               i, j, k, grid, auxiliary, fields, relation, Δt)
        else
            salinity_coupling_flux(i, j, k, grid, C, S)
        end
        left_flux += I[i, j, k]

        right_flux = if k == Nz
            column_top_boundary_energy_flux(boundary_conditions.top,
                                            i, j, k, grid, auxiliary, fields, relation, Δt)
        else
            salinity_coupling_flux(i, j, k+1, grid, C, S)
        end
        right_flux += I[i, j, k+1]

        flux_tendency = (right_flux - left_flux) / current_column_cell_thickness(i, j, k, grid)

        # Conservative moving-grid balance divided by the current layer thickness:
        # Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ).
        moving_flux_difference =
            moving_face_displacement_flux(i, j, k+1, grid, fields.internal_energy,
                                          boundary_conditions.bottom,
                                          boundary_conditions.top) -
            moving_face_displacement_flux(i, j, k, grid, fields.internal_energy,
                                          boundary_conditions.bottom,
                                          boundary_conditions.top)
        moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
        metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] =
            metric_ratio * fields.internal_energy[i, j, k] + moving_tendency + Δt * flux_tendency
    end
end

function _assemble_column_energy_system_cpu!(auxiliary, fields, grid, boundary_conditions, relation, Δt)
    Nz = size(grid, 3)

    for k in 1:Nz, j in 1:size(grid, 2), i in 1:size(grid, 1)
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
                column_bottom_boundary_energy_factor(boundary_conditions.bottom,
                                                     i, j, k, grid, auxiliary, fields, Δt)
            else
                zero(Δt)
            end

            top_factor = if k == Nz
                column_top_boundary_energy_factor(boundary_conditions.top,
                                                  i, j, k, grid, auxiliary, fields, Δt)
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
                column_bottom_boundary_energy_flux(boundary_conditions.bottom,
                                                   i, j, k, grid, auxiliary, fields, relation, Δt)
            else
                salinity_coupling_flux(i, j, k, grid, C, S)
            end
            left_flux += I[i, j, k]

            right_flux = if k == Nz
                column_top_boundary_energy_flux(boundary_conditions.top,
                                                i, j, k, grid, auxiliary, fields, relation, Δt)
            else
                salinity_coupling_flux(i, j, k+1, grid, C, S)
            end
            right_flux += I[i, j, k+1]

            flux_tendency = (right_flux - left_flux) / current_column_cell_thickness(i, j, k, grid)

            # Conservative moving-grid balance divided by the current layer thickness:
            # Jⁿ⁺¹ Δr Eⁿ⁺¹ = Jⁿ Δr Eⁿ + δzᵤ Eᵘᵖᵤ - δzₗ Eᵘᵖₗ + Δt (Fᵤ - Fₗ).
            moving_flux_difference =
                moving_face_displacement_flux(i, j, k+1, grid, fields.internal_energy,
                                              boundary_conditions.bottom,
                                              boundary_conditions.top) -
                moving_face_displacement_flux(i, j, k, grid, fields.internal_energy,
                                              boundary_conditions.bottom,
                                              boundary_conditions.top)
            moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
            metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
            auxiliary.energy_rhs[i, j, k] =
                metric_ratio * fields.internal_energy[i, j, k] + moving_tendency + Δt * flux_tendency
        end
    end

    return nothing
end

"""
    assemble_column_energy_system!(thermodynamics, dt)

Assemble the tridiagonal backward-Euler system for one internal-energy step.
"""
function assemble_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    grid = thermodynamics.fields.internal_energy.grid
    arch = architecture(grid)

    compute_column_shortwave_flux!(thermodynamics)

    if on_cpu(grid)
        _assemble_column_energy_system_cpu!(thermodynamics.auxiliary,
                                            thermodynamics.fields,
                                            grid,
                                            thermodynamics.boundary_conditions,
                                            thermodynamics.relation,
                                            Δt)
    else
        launch!(arch, grid, :xyz,
                _assemble_column_energy_system!,
                thermodynamics.auxiliary,
                thermodynamics.fields,
                grid,
                thermodynamics.boundary_conditions,
                thermodynamics.relation,
                Δt)
    end

    return nothing
end

"""
    solve_column_energy_system!(thermodynamics)

Solve the assembled column energy system into the internal-energy field.
"""
function solve_column_tridiagonal_system_cpu!(solution, solver, rhs)
    grid = solver.grid
    a = solver.a
    b = solver.b
    c = solver.c
    scratch = solver.t
    Nz = size(grid, 3)
    FT = eltype(grid)

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            β = b[i, j, 1]
            solution[i, j, 1] = rhs[i, j, 1] / β

            for k in 2:Nz
                cᵏ⁻¹ = c[i, j, k-1]
                bᵏ = b[i, j, k]
                aᵏ⁻¹ = a[i, j, k-1]

                scratch[i, j, k] = cᵏ⁻¹ / β
                β = bᵏ - aᵏ⁻¹ * scratch[i, j, k]

                ϕ★ = (rhs[i, j, k] - aᵏ⁻¹ * solution[i, j, k-1]) / β
                definitely_diagonally_dominant = abs(β) > 10 * eps(FT)
                solution[i, j, k] = ifelse(definitely_diagonally_dominant,
                                           ϕ★,
                                           solution[i, j, k])
            end

            for k in Nz-1:-1:1
                solution[i, j, k] -= scratch[i, j, k+1] * solution[i, j, k+1]
            end
        end
    end

    return nothing
end

function solve_column_energy_system!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid

    if on_cpu(grid)
        solve_column_tridiagonal_system_cpu!(thermodynamics.fields.internal_energy,
                                             thermodynamics.solvers.energy_solver,
                                             thermodynamics.auxiliary.energy_rhs)
    else
        solve!(thermodynamics.fields.internal_energy,
               thermodynamics.solvers.energy_solver,
               thermodynamics.auxiliary.energy_rhs)
    end

    return nothing
end

column_prescribed_temperature(boundary::PrescribedTemperature{<:Number}, i, j) =
    boundary.temperature

column_prescribed_temperature(boundary::PrescribedTemperature, i, j) =
    @inbounds boundary.temperature[i, j]

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
    icepack_temperature_matrix_step!(thermodynamics, dt; max_iterations=100, tolerance=sqrt(eps(eltype(grid))))

Advance a fixed-salinity column with the Icepack BL99 temperature-matrix
linearization. This validation utility is source-traceable to Icepack
`temperature_changes`: conductance is held at the start-of-step temperature,
and the brine-pocket heat capacity uses the secant form between the initial and
latest predicted temperatures. The supported boundary configuration is a
prescribed top energy flux and prescribed bottom temperature.
"""
function icepack_temperature_matrix_step!(thermodynamics::ColumnEnergyThermodynamics,
                                          Δt;
                                          max_iterations = 100,
                                          tolerance = sqrt(eps(eltype(thermodynamics.fields.internal_energy.grid))))
    grid = thermodynamics.fields.internal_energy.grid
    on_cpu(grid) || error("icepack_temperature_matrix_step! currently supports CPU grids only.")
    thermodynamics.relation isa FixedSalinityBrinePocketEnergyRelation ||
        error("icepack_temperature_matrix_step! requires FixedSalinityBrinePocketEnergyRelation.")
    thermodynamics.boundary_conditions.top isa PrescribedEnergyFlux ||
        error("icepack_temperature_matrix_step! requires PrescribedEnergyFlux at the top boundary.")
    thermodynamics.boundary_conditions.bottom isa PrescribedTemperature ||
        error("icepack_temperature_matrix_step! requires PrescribedTemperature at the bottom boundary.")

    fields = thermodynamics.fields
    auxiliary = thermodynamics.auxiliary
    relation = thermodynamics.relation
    conductivity = thermal_conductivity(thermodynamics.energy_transport)
    top_flux = column_boundary_energy_flux(thermodynamics.boundary_conditions.top)
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
            Δz = current_column_cell_thickness(i, j, 1, grid)
            Tbot = column_prescribed_temperature(thermodynamics.boundary_conditions.bottom, i, j)

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
        moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] = metric_ratio * fields.bulk_salinity[i, j, k] + moving_tendency
    end
end

function _assemble_column_salinity_system_cpu!(auxiliary, fields, grid, Δt)
    Nz = size(grid, 3)

    for k in 1:Nz, j in 1:size(grid, 2), i in 1:size(grid, 1)
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
            moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
            auxiliary.energy_rhs[i, j, k] = metric_ratio * fields.bulk_salinity[i, j, k] + moving_tendency
        end
    end

    return nothing
end

"""
    assemble_column_salinity_system!(thermodynamics, dt)

Assemble the tridiagonal backward-Euler system for one closed-boundary bulk
salinity diffusion step.
"""
function assemble_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    grid = thermodynamics.fields.bulk_salinity.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _assemble_column_salinity_system_cpu!(thermodynamics.auxiliary,
                                              thermodynamics.fields,
                                              grid,
                                              Δt)
    else
        launch!(arch, grid, :xyz,
                _assemble_column_salinity_system!,
                thermodynamics.auxiliary,
                thermodynamics.fields,
                grid,
                Δt)
    end

    return nothing
end

"""
    solve_column_salinity_system!(thermodynamics)

Solve the assembled scalar salinity system into the bulk-salinity field.
"""
function solve_column_salinity_system!(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.bulk_salinity.grid

    if on_cpu(grid)
        solve_column_tridiagonal_system_cpu!(thermodynamics.fields.bulk_salinity,
                                             thermodynamics.solvers.salinity_solver,
                                             thermodynamics.auxiliary.energy_rhs)
    else
        solve!(thermodynamics.fields.bulk_salinity,
               thermodynamics.solvers.salinity_solver,
               thermodynamics.auxiliary.energy_rhs)
    end

    return nothing
end

@kernel function _compute_column_metric_change!(scratch, field, grid)
    i, j, k = @index(Global, NTuple)
    metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
    moving_flux_difference =
        moving_face_displacement_flux(i, j, k+1, grid, field) -
        moving_face_displacement_flux(i, j, k, grid, field)
    moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
    @inbounds scratch[i, j, k] = metric_ratio * field[i, j, k] + moving_tendency
end

@kernel function _copy_column_metric_change!(field, scratch)
    i, j, k = @index(Global, NTuple)
    @inbounds field[i, j, k] = scratch[i, j, k]
end

function _apply_column_metric_change_cpu!(field, scratch)
    grid = field.grid

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            metric_ratio = previous_to_current_column_metric_ratio(i, j, k, grid)
            moving_flux_difference =
                moving_face_displacement_flux(i, j, k+1, grid, field) -
                moving_face_displacement_flux(i, j, k, grid, field)
            moving_tendency = moving_flux_difference / current_column_cell_thickness(i, j, k, grid)
            scratch[i, j, k] = metric_ratio * field[i, j, k] + moving_tendency
        end
    end

    for k in 1:size(grid, 3), j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds field[i, j, k] = scratch[i, j, k]
    end

    return nothing
end

function apply_column_metric_change!(field, scratch)
    grid = field.grid
    arch = architecture(grid)

    if on_cpu(grid)
        _apply_column_metric_change_cpu!(field, scratch)
    else
        launch!(arch, grid, :xyz, _compute_column_metric_change!, scratch, field, grid)
        launch!(arch, grid, :xyz, _copy_column_metric_change!, field, scratch)
    end

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
        ice_thickness[i, j, 1] += column_stefan_thickness_change(phase_transitions,
                                                                 ρi,
                                                                 δQ,
                                                                 Δt)
    end
end

function _column_stefan_thickness_update_cpu!(ice_thickness,
                                              phase_transitions,
                                              sea_ice_density,
                                              residual_energy_flux,
                                              Δt)
    grid = ice_thickness.grid

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        @inbounds begin
            ρi = column_stefan_parameter(sea_ice_density, i, j, 1)
            δQ = column_stefan_parameter(residual_energy_flux, i, j, 1)
            ice_thickness[i, j, 1] += column_stefan_thickness_change(phase_transitions,
                                                                     ρi,
                                                                     δQ,
                                                                     Δt)
        end
    end

    return nothing
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

    if on_cpu(grid)
        _column_stefan_thickness_update_cpu!(ice_thickness,
                                             phase_transitions,
                                             sea_ice_density,
                                             residual_energy_flux,
                                             Δt)
    else
        launch!(arch, grid, :xy,
                _column_stefan_thickness_update!,
                ice_thickness,
                phase_transitions,
                sea_ice_density,
                residual_energy_flux,
                Δt)
    end

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

"""
    column_energy_time_step!(thermodynamics, dt)

Advance the column thermodynamics by one semi-implicit energy step,
followed by the scalar salinity step when salinity is prognostic. On a
`MutableVerticalDiscretization`, the energy and salinity updates include the
previous-to-current metric ratio and explicit swept-face tracer fluxes.
"""
function column_energy_time_step!(thermodynamics::ColumnEnergyThermodynamics, Δt)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    compute_column_transport_coefficients!(thermodynamics)
    assemble_column_energy_system!(thermodynamics, Δt)
    solve_column_energy_system!(thermodynamics)
    compute_column_thermodynamic_diagnostics!(thermodynamics)
    column_salinity_time_step!(thermodynamics, Δt)

    return nothing
end

function thermodynamic_time_step!(model,
                                  thermodynamics::ColumnEnergyThermodynamics,
                                  ::Nothing,
                                  Δt)
    column_energy_time_step!(thermodynamics, Δt)
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
        total_energy += @inbounds E[i, j, k] * current_column_cell_thickness(i, j, k, grid)
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
        total_salinity += @inbounds S[i, j, k] * current_column_cell_thickness(i, j, k, grid)
    end

    return total_salinity
end

function column_boundary_energy_flux_difference(thermodynamics::ColumnEnergyThermodynamics)
    grid = thermodynamics.fields.internal_energy.grid
    boundary_conditions = thermodynamics.boundary_conditions
    difference = zero(eltype(grid))
    boundary_difference = (column_boundary_energy_flux(boundary_conditions.top) -
                           column_boundary_energy_flux(boundary_conditions.bottom))

    for j in 1:size(grid, 2), i in 1:size(grid, 1)
        difference += boundary_difference
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
    column_energy_budget(thermodynamics, initial_energy, dt; final_energy = column_integrated_energy(thermodynamics),
                         surface_stefan_residual_flux = 0)

Return a named tuple summarizing the column-integrated energy budget over one
step with duration `dt`.
"""
function column_energy_budget(thermodynamics::ColumnEnergyThermodynamics,
                              initial_energy,
                              Δt;
                              final_energy = column_integrated_energy(thermodynamics),
                              surface_stefan_residual_flux = 0)
    boundary_flux_change = Δt * column_boundary_energy_flux_difference(thermodynamics)
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

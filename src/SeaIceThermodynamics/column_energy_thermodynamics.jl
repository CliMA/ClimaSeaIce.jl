import Oceananigans: prognostic_state, restore_prognostic_state!
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: CPU, architecture
using Oceananigans.Fields: AbstractField
using Oceananigans.Grids: ZDirection, znode
using Oceananigans.Operators: Δzᵃᵃᶜ, Δzᵃᵃᶠ
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
    ConductiveTemperatureTransport(convert(FT, conductivity))

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
    ConductiveAndDiffusiveEnergyTransport(convert(FT, conductivity),
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
    MeltingLimitedSurfaceFlux([FT=Oceananigans.defaults.FloatType; flux=0])

Top energy boundary condition with an imposed surface flux that is capped by
the complete-melt energy of the top cell. The applied flux warms the fixed-grid
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

Top and bottom energy boundary conditions for the fixed-grid column energy
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

function initialize_salinity!(fields, salinity_closure)
    profile = salinity_profile(salinity_closure)
    isnothing(profile) || set!(fields.bulk_salinity, profile)
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

    isnothing(internal_energy) || set!(thermodynamic_fields.internal_energy, internal_energy)
    isnothing(bulk_salinity) || set!(thermodynamic_fields.bulk_salinity, bulk_salinity)
    isnothing(temperature) || set!(thermodynamic_fields.temperature, temperature)
    isnothing(liquid_fraction) || set!(thermodynamic_fields.liquid_fraction, liquid_fraction)
    isnothing(brine_salinity) || set!(thermodynamic_fields.brine_salinity, brine_salinity)

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
##### Diagnostics and fixed-grid energy step
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

@kernel function _compute_column_transport_coefficients!(auxiliary, fields, energy_transport)
    i, j, k = @index(Global, NTuple)
    Nz = size(fields.internal_energy.grid, 3)

    kappa_T = thermal_conductivity(energy_transport)
    kappa_E = energy_diffusivity(energy_transport)

    @inbounds begin
        if k == 1
            auxiliary.thermal_conductivity[i, j, 1] = zero(kappa_T)
            auxiliary.energy_diffusivity[i, j, 1] = zero(kappa_E)
            auxiliary.effective_energy_diffusivity[i, j, 1] = zero(kappa_T)
            auxiliary.salinity_coupling_diffusivity[i, j, 1] = zero(kappa_T)
        end

        if k > 1
            TE = (fields.temperature_energy_derivative[i, j, k-1] +
                  fields.temperature_energy_derivative[i, j, k]) / 2

            TS = (fields.temperature_salinity_derivative[i, j, k-1] +
                  fields.temperature_salinity_derivative[i, j, k]) / 2

            auxiliary.thermal_conductivity[i, j, k] = kappa_T
            auxiliary.energy_diffusivity[i, j, k] = kappa_E
            auxiliary.effective_energy_diffusivity[i, j, k] = kappa_T * TE + kappa_E
            auxiliary.salinity_coupling_diffusivity[i, j, k] = kappa_T * TS
        end

        if k == Nz
            kp1 = Nz + 1
            auxiliary.thermal_conductivity[i, j, kp1] = zero(kappa_T)
            auxiliary.energy_diffusivity[i, j, kp1] = zero(kappa_E)
            auxiliary.effective_energy_diffusivity[i, j, kp1] = zero(kappa_T)
            auxiliary.salinity_coupling_diffusivity[i, j, kp1] = zero(kappa_T)
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
            auxiliary.thermal_conductivity[i, j, 1] = zero(kappa_T)
            auxiliary.energy_diffusivity[i, j, 1] = zero(kappa_E)
            auxiliary.effective_energy_diffusivity[i, j, 1] = zero(kappa_T)
            auxiliary.salinity_coupling_diffusivity[i, j, 1] = zero(kappa_T)

            for k in 2:Nz
                TE = (fields.temperature_energy_derivative[i, j, k-1] +
                      fields.temperature_energy_derivative[i, j, k]) / 2

                TS = (fields.temperature_salinity_derivative[i, j, k-1] +
                      fields.temperature_salinity_derivative[i, j, k]) / 2

                auxiliary.thermal_conductivity[i, j, k] = kappa_T
                auxiliary.energy_diffusivity[i, j, k] = kappa_E
                auxiliary.effective_energy_diffusivity[i, j, k] = kappa_T * TE + kappa_E
                auxiliary.salinity_coupling_diffusivity[i, j, k] = kappa_T * TS
            end

            auxiliary.thermal_conductivity[i, j, Nz+1] = zero(kappa_T)
            auxiliary.energy_diffusivity[i, j, Nz+1] = zero(kappa_E)
            auxiliary.effective_energy_diffusivity[i, j, Nz+1] = zero(kappa_T)
            auxiliary.salinity_coupling_diffusivity[i, j, Nz+1] = zero(kappa_T)
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
        Δz = Δzᵃᵃᶜ(i, j, Nz, grid)

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
            Δz = Δzᵃᵃᶜ(i, j, Nz, grid)

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
            Δz = Δzᵃᵃᶜ(i, j, Nz, grid)

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
           (Δzᵃᵃᶜ(i, j, kc, grid) * Δzᵃᵃᶠ(i, j, kf, grid))
end

@inline function salinity_coupling_flux(i, j, k, grid, C, S)
    return @inbounds(C[i, j, k] * (S[i, j, k] - S[i, j, k-1]) / Δzᵃᵃᶠ(i, j, k, grid))
end

@inline column_boundary_energy_flux(::InsulatingBoundary) = 0
@inline column_boundary_energy_flux(boundary::PrescribedEnergyFlux) = boundary.flux
@inline column_boundary_energy_flux(boundary::MeltingLimitedSurfaceFlux) = boundary.flux

@inline column_requested_surface_energy_flux(boundary) = column_boundary_energy_flux(boundary)

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
    Δz = Δzᵃᵃᶜ(i, j, k, grid)
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

        left_flux = if k == 1
            column_boundary_energy_flux(boundary_conditions.bottom, i, j, k, grid, fields, relation, Δt)
        else
            salinity_coupling_flux(i, j, k, grid, C, S)
        end
        left_flux += I[i, j, k]

        right_flux = if k == Nz
            column_boundary_energy_flux(boundary_conditions.top, i, j, k, grid, fields, relation, Δt)
        else
            salinity_coupling_flux(i, j, k+1, grid, C, S)
        end
        right_flux += I[i, j, k+1]

        flux_tendency = (right_flux - left_flux) / Δzᵃᵃᶜ(i, j, k, grid)
        auxiliary.energy_rhs[i, j, k] = fields.internal_energy[i, j, k] + Δt * flux_tendency
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

            left_flux = if k == 1
                column_boundary_energy_flux(boundary_conditions.bottom, i, j, k, grid, fields, relation, Δt)
            else
                salinity_coupling_flux(i, j, k, grid, C, S)
            end
            left_flux += I[i, j, k]

            right_flux = if k == Nz
                column_boundary_energy_flux(boundary_conditions.top, i, j, k, grid, fields, relation, Δt)
            else
                salinity_coupling_flux(i, j, k+1, grid, C, S)
            end
            right_flux += I[i, j, k+1]

            flux_tendency = (right_flux - left_flux) / Δzᵃᵃᶜ(i, j, k, grid)
            auxiliary.energy_rhs[i, j, k] = fields.internal_energy[i, j, k] + Δt * flux_tendency
        end
    end

    return nothing
end

"""
    assemble_column_energy_system!(thermodynamics, dt)

Assemble the tridiagonal backward-Euler system for one fixed-grid internal
energy step.
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

        auxiliary.energy_rhs[i, j, k] = fields.bulk_salinity[i, j, k]
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

            auxiliary.energy_rhs[i, j, k] = fields.bulk_salinity[i, j, k]
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

"""
    column_salinity_time_step!(thermodynamics, dt)

Advance prognostic bulk salinity by one scalar fixed-grid diffusion step. The
step is a no-op for prescribed salinity and for [`NoSalinityTransport`](@ref).
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

"""
    column_energy_time_step!(thermodynamics, dt)

Advance the column thermodynamics by one fixed-grid semi-implicit energy step,
followed by the scalar salinity step when salinity is prognostic.
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
        total_energy += @inbounds E[i, j, k] * Δzᵃᵃᶜ(i, j, k, grid)
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
        total_salinity += @inbounds S[i, j, k] * Δzᵃᵃᶜ(i, j, k, grid)
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

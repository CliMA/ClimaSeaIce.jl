using ClimaSeaIce.Rheologies

struct SeaIceMomentumEquation{S, C, R, A}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
    solver :: S
end

struct ExplicitSolver end

function SeaIceMomentumEquation(grid; 
                                coriolis = nothing,
                                rheology = ElastoViscoPlasticRheology(eltype(grid)),
                                auxiliary_fields = NamedTuple(),
                                solver = ExplicitSolver())

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))

    return SeaIceMomentumEquation(coriolis, rheology, auxiliary_fields, solver)
end

fields(mom::SeaIceMomentumEquation) = mom.auxiliary_fields
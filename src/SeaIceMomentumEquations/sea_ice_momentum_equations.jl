struct SeaIceMomentumEquation{S, C, R, A}
    coriolis :: C
    rheology :: R
    auxiliary_fields :: A
    solver :: S
end

struct ExplicitSolver end

function SeaIceMomentumEquation(grid; 
                                  coriolis = nothing,
                                  rheology = nothing,
                                  auxiliary_fields = NamedTuple(),
                                  solver = ExplicitSolver())

    auxiliary_fields = merge(auxiliary_fields, required_auxiliary_fields(rheology, grid))

    return SeaIceMomentumEquation(coriolis, rheology, auxiliary_fields, solver)
end

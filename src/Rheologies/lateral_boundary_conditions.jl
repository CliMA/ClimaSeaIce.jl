using Oceananigans.BoundaryConditions: BoundaryCondition, Value
using Oceananigans.Grids: Center, Face
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, ImmersedBoundaryCondition,
                                       immersed_inactive_node, immersed_peripheral_node

struct FreeSlip end
struct NoSlip end

const VBC = BoundaryCondition{<:Value}

# `ImmersedBoundaryCondition` orders its sides (west, east, south, north, bottom, top)
const NoSlipU = ImmersedBoundaryCondition{<:Any, <:Any, <:VBC, <:VBC}
const NoSlipV = ImmersedBoundaryCondition{<:VBC, <:VBC}

@inline u_lateral_boundary_condition(immersed_boundary_condition) = FreeSlip()
@inline u_lateral_boundary_condition(::NoSlipU) = NoSlip()

@inline v_lateral_boundary_condition(immersed_boundary_condition) = FreeSlip()
@inline v_lateral_boundary_condition(::NoSlipV) = NoSlip()

# The shear strain rate mixes ∂u/∂y and ∂v/∂x, so it is no-slip if either component is.
@inline combine_slip(::FreeSlip, ::FreeSlip) = FreeSlip()
@inline combine_slip(::NoSlip,   ::FreeSlip) = NoSlip()
@inline combine_slip(::FreeSlip, ::NoSlip)   = NoSlip()
@inline combine_slip(::NoSlip,   ::NoSlip)   = NoSlip()

@inline lateral_boundary_condition(u_immersed_bc, v_immersed_bc) = 
    combine_slip(u_lateral_boundary_condition(u_immersed_bc),
                 v_lateral_boundary_condition(v_immersed_bc))

@inline strain_rate_slip_factor(::FreeSlip) = 1
@inline strain_rate_slip_factor(::NoSlip) = 2

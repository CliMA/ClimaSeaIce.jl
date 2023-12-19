using Base: @propagate_inbounds

import Adapt
import Base: getindex

"""
    DifferenceOfArrays{F}

`DifferenceOfArrays` objects hold 2 arrays/fields and return their difference when indexed.
"""
struct DifferenceOfArrays{F}
    arrays :: F
end

DifferenceOfArrays(array1, array2) = new((array1, array2))

@propagate_inbounds getindex(s::DifferenceOfArrays, i...) = 
    getindex(s.arrays[1], i...) - getindex(s.arrays[2], i...)

Adapt.adapt_structure(to, s::DifferenceOfArrays) = DifferenceOfArrays(Adapt.adapt(to, s.arrays))
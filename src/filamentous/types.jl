
# Type definitions such that we can have Articulation{FilamentProperties},
# allowing torque_elastic to be applied to an articulation.

abstract type ElasticProperties end;

mutable struct LinearElasticProperties <: ElasticProperties
    spatial_rigidity::Matrix{<:Real} #3x3 matrix of elastic rigidity.
end

mutable struct FilamentProperties
    elastics::ElasticProperties
end

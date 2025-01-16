

abstract type ElasticProperties end;

mutable struct LinearElasticProperties <: ElasticProperties
    spatial_rigidity::Matrix{<:Real}
end

mutable struct FilamentProperties
    elastics::ElasticProperties
end
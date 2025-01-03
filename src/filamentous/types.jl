

abstract type ElasticProperties end;

mutable struct LinearElasticProperties <: ElasticProperties
    spatial_rigidity::Vector{<:Real}
end
AxisymmetricElasticProperties(spatial_rigidity::Real) = LinearElasticProperties([spatial_rigidity, spatial_rigidity, 0, 0, 0, 0]);

mutable struct FilamentProperties
    elastics::ElasticProperties
end
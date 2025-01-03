

torque_elastic() = InternalTorque((a::Articulation, i::ArticulationHarness, t::Real) -> elastic_response(a.joint_type, a.properties.elastics, i.q));

function elastic_response(j::RotaryJoint, properties::LinearElasticProperties, theta::Vector{<:Real})
    if j.axis == :x
        return -properties.spatial_rigidity[1] * theta;
    end
end

# The driving torque for elasticity. Requires elastic_response to be defined for the joint type
# and requires that every articulation.properties has a .elastics.
torque_elastic() = InternalTorque((a::Articulation, i::ArticulationHarness, t::Real) -> elastic_response(a.joint_type, a.properties.elastics, i));


function elastic_response(j::FixedJoint, properties::LinearElasticProperties, i::ArticulationHarness)
    return zeros(0);
end

function elastic_response(j::RotaryJoint, properties::LinearElasticProperties, i::ArticulationHarness)
    if j.axis == :x
        return -properties.spatial_rigidity[1,1] * i.q;
    elseif j.axis == :y
        return -properties.spatial_rigidity[2,2] * i.q;
    else
        return -properties.spatial_rigidity[3,3] * i.q;
    end
end

function elastic_response(j::EulerXYJoint, properties::LinearElasticProperties, i::ArticulationHarness)
    EI_l = properties.spatial_rigidity;
    return -[EI_l[1,1],EI_l[2,2]] .* i.q;
end

function elastic_response(j::BendingJoint, properties::LinearElasticProperties, i::ArticulationHarness)
    
    # e_z = [0,0,1];
    qx = i.q[1]; qy = i.q[2];
    sx = sin(qx); sy = sin(qy);
    cx = cos(qx); cy = cos(qy);

    # e_d = [sx*cy, -sy, cx*cy];

    torque_axis = [sy, sx*cy];
    # = [ -e_d[2], e_d[1], 0];
    # = cross(undeformed_vector, deformed_vector);

    EI = properties.spatial_rigidity[1:2,1:2];

    ed_z = cx*cy # = e_d[3];
    acos_edz = acos(ed_z);
    gamma = acos_edz / sin(acos_edz);
    torque = (EI * torque_axis) * ifelse(ed_z==1, 1, gamma);
    bending = i.S[1:2,1:2]' * torque;
    
    return bending;
    
end


# Return a vector of stiffnesses of each DOF.
function get_stiffnesses(body::AbstractArticulatedBody)
    stiffnesses = zeros(dof(body));
    forward_recurse_iterative!(body, get_stiffness, stiffnesses);
    return stiffnesses;
end
function get_stiffness(a::Articulation, stiffnesses::Vector{<:Real})
    individual_harness = ArticulationHarness(a, Float64);
    XJ, vJ, cJ, S = j_calc(a.joint_type, individual_harness);
    stiffnesses[a.state_indices] = diag(S[1:3,:]' * a.properties.elastics.spatial_rigidity * S[1:3,:]);
end

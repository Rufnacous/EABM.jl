
# See Chapter 4.4 of Rigid Body Dynamics Algorithms by Roy Featherstone (2008).
# Look for "jcalc" or "joint calculation routine".
function j_calc(j::JointType, s::ArticulationHarness)
    S,Ṡ,E,p = get_transformations(j, s.q, s.q_dt);

    numerictype = typeof(E).parameters[1];
    # Copying the numerictype allows this to be auto-differentiable by ForwardDiff.
    E_diag = zeros(numerictype, 6, 6);
    E_diag[1:3,1:3] = E; E_diag[4:6,4:6] = E;
    p_diag = xlt(p);

    XJ = E_diag * p_diag

    vJ = Vector{numerictype}(S * s.q_dt); 
    cJ = Ṡ * s.q_dt;

    return XJ, vJ, cJ, S
end

# Just to be use for the articulation_zero of a body.
struct ArticulationZeroJoint <: JointType end
dof(j::ArticulationZeroJoint) = 0;



# A joint which rigidly connects an articulation to its parent. i.e. with no articulation at all.
struct FixedJoint <: JointType
end
dof(j::FixedJoint) = 0;
function get_transformations(j::FixedJoint, q::Vector{<: Real}, qdt::Vector{<: Real})
    return (
            zeros(6,0),
            zeros(6,0),
            [1 0 0; 0 1 0; 0 0 1],
            [0,0,0]
        )
end


# 1DOF rotational joints. RotaryJoint(:x), RotaryJoint(:y), RotaryJoint(:z).
struct RotaryJoint <: JointType
    axis::Symbol
end
dof(j::RotaryJoint) = 1;
function get_transformations(j::RotaryJoint, q::Vector{<: Real}, qdt::Vector{<: Real})
    if j.axis == :x
        return (
            Matrix([1 0 0 0 0 0]'),
            Matrix([0 0 0 0 0 0]'),
            rotate_x(q[1]),
            [0,0,0]
        )
    elseif j.axis ==:y
        return (
            Matrix([0 1 0 0 0 0]'),
            Matrix([0 0 0 0 0 0]'),
            rotate_y(q[1]),
            [0,0,0]
        )
    else
        return (
            Matrix([0 0 1 0 0 0]'),
            Matrix([0 0 0 0 0 0]'),
            rotate_z(q[1]),
            [0,0,0]
        )
    end
end

# 2DOF bending joints akin to Euler Rotations. Suffer from small-angle limitations.
struct EulerXYJoint <: JointType end
dof(j::EulerXYJoint) = 2;
function get_transformations(j::EulerXYJoint, q::Vector{<: Real}, qdt::Vector{<: Real})
    c1 = cos(0); s1 = sin(0); c2 = cos(q[1]); s2 = sin(q[1]); c3 = cos(q[2]); s3 = sin(q[2]);
    E = [
        (c2)    (0)   (-s2)   ;
        (s2*s3) (c3)  (c2*s3) ;
        (s2*c3) (-s3) (c2*c3) ];
    S = [
        0     1;
        c3    0;
        -s3   0;
        0     0;
        0     0;
        0     0
    ]
    Ṡ = [
        0               0;
        -s3*qdt[2]       0;
        -c3*qdt[2]       0;
        0               0;
        0               0;
        0               0
    ]
    p = [0,0,0];
    return (S, Ṡ, E, p);
end



# Better 2DOF bending joints.
struct BendingJoint <: JointType end
dof(j::BendingJoint) = 2;
function get_transformations(j::BendingJoint, q::Vector{<: Real}, qdt::Vector{<: Real})

    cy = cos(q[1]); sy = sin(q[1]); cz = cos(q[2]); sz = sin(q[2]);
    qdtz = qdt[2];

    E = get_rotation_b_is_ez([sy*cz, -sz, cy*cz]);

    S = [ 
        0     1;
        cz    0;
        0     0;
        0     0;
        0     0;
        0     0
    ]
    Ṡ = [
        0               0;
        -sz*qdtz        0;
        0               0;
        0               0;
        0               0;
        0               0
    ]
    p = [0,0,0];
    return (S, Ṡ, E, p);

end

function j_calc(j::JointType, s::ArticulationHarness)
    S,Ṡ,E,p = get_transformations(j, s.q, s.q_dt);

    numerictype = typeof(E).parameters[1];
    E_diag = zeros(numerictype, 6, 6);
    E_diag[1:3,1:3] = E; E_diag[4:6,4:6] = E;
    p_diag = xlt(p);

    XJ = E_diag * p_diag

    vJ = Vector{numerictype}(S * s.q_dt); 
    cJ = Ṡ * s.q_dt;

    return XJ, vJ, cJ, S
end

struct ArticulationZeroJoint <: JointType end
dof(j::ArticulationZeroJoint) = 0;

struct RotaryJoint <: JointType
    axis::Symbol
end
dof(j::RotaryJoint) = 1;
function get_transformations(j::RotaryJoint, q::Vector{<: Real}, qdt::Vector{<: Real})
    c1 = cos(q[1]); s1 = sin(q[1]);
    return (
        Matrix([1 0 0 0 0 0]'),
        Matrix([0 0 0 0 0 0]'),
        [1 0 0; 0 c1 s1; 0 -s1 c1],
        [0,0,0]
    )
    # if j.axis == :x
    #     return (
    #         Matrix([1 0 0 0 0 0]'),
    #         Matrix([0 0 0 0 0 0]'),    # Ṡ
    #         [ 1   0   0;
    #           0  c1  s1;
    #           0 -s1  c1],              # E
    #         [0,0,0]                    # p
    #     )
end
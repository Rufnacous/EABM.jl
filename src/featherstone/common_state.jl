
function get_acceleration(body::AbstractArticulatedBody, harness::StateHarness)
    # returns the accelerations of an articulated body as a vector.
    accel = zeros(typeof(harness.articulation_states[end].q_d²t).parameters[1], dof(body));
    forward_recurse_iterative!(body, get_acceleration, accel, harness);
    return accel;
end
function get_acceleration(a::Articulation, accel::Vector{<:Real}, harness::StateHarness)
    # inner loop function for the above.
    accel[a.state_indices] = harness[a].q_d²t;
end

function get_acceleration!(body::AbstractArticulatedBody, harness::StateHarness, accel::Vector{<:Real})
    # returns the accelerations of an articulated body as a vector.
    forward_recurse_iterative!(body, get_acceleration!, accel, harness);
end
function get_acceleration!(a::Articulation, accel::Vector{<:Real}, harness::StateHarness)
    # inner loop function for the above.
    accel[a.state_indices] = harness[a].q_d²t;
end





function get_position(body::AbstractArticulatedBody, state::Vector{<:Real}; with_retrace::Bool=false)
    # returns the cartesian positions of articulations in a 3xN matrix. with_retrace is useful for line plotting.
    if with_retrace
        pos = [];
        harness = StateHarness(state, body);
        push!(pos, harness.articulation_states[1].p);
        forward_recurse!(body, push_position, pos, harness, retrace=true);
        
        pos_out = zeros(3, 2+length(pos));
        pos_out[:,1] = harness.articulation_states[1].p;
        for i in eachindex(pos)
            pos_out[:,i+1] = pos[i];
        end
        pos_out[:,end] = harness.articulation_states[1].p;

        return pos_out
    else
        pos = zeros(3, 1+n_bodies(body));
        harness = StateHarness(state, body);
        pos[:,1] = harness.articulation_states[1].p;
        forward_recurse!(body, get_position, pos, harness);
        return pos;
    end
end
function get_position(a::Articulation, pos::Matrix{<:Real}, harness::StateHarness)
    # inner loop function for the above.
    pos[:, 1+a.body_number] = harness[a].p;
end
function push_position(a::Articulation, pos, harness::StateHarness)
    # also an inner loop function for get_position.
    push!(pos, harness[a].p);
end






function Base.zeros(articulated_body::AbstractArticulatedBody; no_derivatives::Bool=false)
    # useful for creating initial conditions.
    return zeros(ifelse(no_derivatives,1,2) * dof(articulated_body));
end





function set_state!(articulated_body::AbstractArticulatedBody, harness::StateHarness, generalized_state::Vector{<: Real})
    # loads a state vector into a state harness, for use by the ABA.
    forward_recurse_iterative!(articulated_body, set_state!, articulated_body, harness, generalized_state, length(generalized_state) == 2dof(articulated_body));
end
function set_state!(a::Articulation, articulated_body::AbstractArticulatedBody, harness::StateHarness, generalized_state::Vector{<: Real}, includes_velocities::Bool)
    q = generalized_state[a.state_indices]; 
    if includes_velocities
        qdt = generalized_state[dof(articulated_body) .+ a.state_indices]; 
        set_state!(a, harness, q, qdt);
    else
        set_state!(a, harness, q, zeros(length(a.state_indices)));
    end
end

function set_state!(a::Articulation, harness::StateHarness, q::AbstractVector{<: Real}, q_dt::AbstractVector{<: Real})
    # Processes the raw states of q and q_dt into the state harness, for velocities, transforms, and cartesian positions.
    i = harness[a];
    λi = harness[λ(a)];

    i.q .= q;
    i.q_dt .= q_dt;

    XJ = i.λX; #@view i.λX[:,:];
    vJ = @view i.λX⁻ᵀ[1,:];
    cJ = @view i.λX⁻ᵀ[2,:];
    j_calc!(a.joint_type, i,
        XJ, vJ, cJ, i.S, i.Ṡ);

    mul!(i.λX⁻ᵀ, a.geometry_transform, a.pregeometry_transform);
    mul!(i.Xλ, XJ, i.λX⁻ᵀ);

    mul!(i.v, i.Xλ, λi.v);
    i.v += vJ;
    spatial_cross!(i.c, i.v, vJ);
    i.c = i.c + cJ;
    
    i.λX = inv(i.Xλ); 

    

    # # OT = i.S * i.q + [0,0,0, (a.child_anchor )...];  Needs reinstating if using prismatic joints.

    G = @view i.λX⁻ᵀ[1:3,1:3];
    G1 = @view λi.X0[4:6,4:6];
    G2 = @view a.pregeometry_transform[4:6,4:6];
    mul!(G, G1', G2');

    mul!(i.p, G, a.pregeometry)
    p_add = @view i.λX⁻ᵀ[3,1:3];
    iX0 = @view i.X0[4:6,4:6];
    mul!(p_add, iX0', a.child_anchor)
    i.p += p_add;
    i.p += λi.p;

    
    xlt_OT = i.X0;
    xlt!(xlt_OT, a.child_anchor);


    
    mul!(i.X0, i.Xλ, λi.X0);
    
    mul!(i.V, i.X0', xlt_OT * i.v)
    i.λX⁻ᵀ .= i.Xλ'; #inv(i.λX)';
end



# Initialise a stateharness with a state vector.
function StateHarness(generalized_state::Vector{<:Real}, body::AbstractArticulatedBody)
    s = StateHarness(typeof(generalized_state).parameters[1], body);
    set_state!(body, s, generalized_state);
    return s;
end

# Initialise a null state harness for a certain type. numeric_type tracking allows use with ForwardDiff.
function StateHarness(numeric_type::DataType, body::AbstractArticulatedBody)
    a_harnesses::Vector{ArticulationHarness} = [];

    push_articulation_harness!(body.articulation_zero, numeric_type, a_harnesses)
    a_harnesses[1].v = [0,0,0,0,0,0];
    a_harnesses[1].p = [0,0,0];
    a_harnesses[1].X0 = Matrix{Float64}(I(6));

    forward_recurse_iterative!(body, push_articulation_harness!, numeric_type, a_harnesses)

    return StateHarness(a_harnesses)
end
function push_articulation_harness!(a::Articulation, numeric_type::DataType, a_harnesses::Vector{ArticulationHarness})
    push!( a_harnesses, ArticulationHarness(a, numeric_type) )
    
end

ArticulationHarness( a::Articulation, numeric_type::DataType ) = 
    ArticulationHarness(
        zeros( numeric_type, dof(a) ),
        zeros( numeric_type, dof(a) ),
        zeros( numeric_type, dof(a) ),

        zeros( numeric_type, 6, 6 ),
        zeros( numeric_type, 6, 6 ),
        zeros( numeric_type, 6, 6 ),

        zeros( numeric_type, 6, 6 ),
        zeros( numeric_type, 6),

        zeros( numeric_type, 6, dof(a) ),
        zeros( numeric_type, dof(a), dof(a) ),
        zeros( numeric_type, dof(a), dof(a) ),
        zeros( numeric_type, dof(a) ),
        
        zeros( numeric_type, 6, dof(a) ),
        zeros( numeric_type, 6, dof(a) ),
        
        zeros( numeric_type, 6 ),
        zeros( numeric_type, 6 ),
        zeros( numeric_type, 6 ),
        
        zeros( numeric_type, 3 ),
        zeros( numeric_type, 6 ),
        zeros( numeric_type, 6, 6 )
    );



# This returns a vector of the u values of each articulation. u is the torque minus the basic load.
function get_torque(body::AbstractArticulatedBody, harness::StateHarness)
    torques = zeros(typeof(harness.articulation_states[end].u).parameters[1], dof(body));
    forward_recurse_iterative!(body, get_torque, torques, harness);
    return torques;
end
function get_torque(a::Articulation, torques::Vector{<:Real}, harness::StateHarness)
    torques[a.state_indices] = harness[a].u;
end
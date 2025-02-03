
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
    if includes_velocities
        set_state!(a, harness, generalized_state[a.state_indices], generalized_state[dof(articulated_body) .+ a.state_indices]);
    else
        set_state!(a, harness, generalized_state[a.state_indices], zeros(length(a.state_indices)));
    end
end

function set_state!(a::Articulation, harness::StateHarness, q::Vector{<: Real}, q_dt::Vector{<: Real})
    # Processes the raw states of q and q_dt into the state harness, for velocities, transforms, and cartesian positions.

    harness[a].q = q;
    harness[a].q_dt = q_dt;

    XJ, vJ, cJ, S = j_calc(a.joint_type, harness[a]);

    harness[a].S = S;
    
    harness[a].Xλ = XJ * a.geometry_transform * a.pregeometry_transform;
    harness[a].λX = inv(harness[a].Xλ); 
    harness[a].λX⁻ᵀ = Matrix(harness[a].Xλ'); #inv(harness[a].λX)';


    harness[a].v = (harness[a].Xλ * harness[λ(a)].v) + vJ;
    harness[a].c = cJ + (harness[a].v ⨱ vJ);

    harness[a].X0 = harness[a].Xλ * harness[λ(a)].X0;
    

    OT = harness[a].S * harness[a].q + [0,0,0, (a.child_anchor )...];  

    harness[a].p = harness[λ(a)].p + 
        (harness[λ(a)].X0[4:6,4:6]' * a.pregeometry_transform[4:6,4:6]' * a.pregeometry) + 
        (harness[a].X0[4:6,4:6]' * OT[4:6]);

    harness[a].V = harness[a].X0' * (xlt(OT[4:6]) * harness[a].v)
end
    



# Initialise a stateharness with a state vector.
function StateHarness(generalized_state::Vector{<:Real}, body::AbstractArticulatedBody)
    s = StateHarness(typeof(generalized_state).parameters[1], body);
    set_state!(body, s, generalized_state);
    return s;
end

# Initialise a null state harness for a certain type. numeric_type tracking allows use with ForwardDiff.
function StateHarness(numeric_type::DataType, body::AbstractArticulatedBody)
    a_harnesses::Vector{ArticulationHarness} = Vector{ArticulationHarness}(undef, 1+n_bodies(body));
    sh = StateHarness(a_harnesses);

    push_articulation_harness!(body.articulation_zero, numeric_type, sh)
    a_harnesses[1].v = [0,0,0,0,0,0];
    a_harnesses[1].p = [0,0,0];
    a_harnesses[1].X0 = Matrix{Float64}(I(6));

    forward_recurse_iterative!(body, push_articulation_harness!, numeric_type, sh)

    return sh
end
function push_articulation_harness!(a::Articulation, numeric_type::DataType, sh::StateHarness)
    sh[a] = ArticulationHarness(a, numeric_type);
    
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
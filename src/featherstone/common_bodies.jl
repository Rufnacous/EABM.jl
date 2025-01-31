
# An example body. No idea why you would use this library to simulate a structure without elasticity.
mutable struct NPendulum <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function NPendulum(n::Integer, seg_lengths::Number, seg_masses::Number, joints::JointType; compound::Symbol=:cylinder, compoundargs...)
    a0 = articulation_zero();
    alast =  a0; anext = nothing;
    next_free_state_index = 1;
    
    for ni = 1:n
        if compound == :cylinder
            anext = CylinderArticulation(ni, next_free_state_index, joints, alast, seg_lengths, compoundargs[:radius], seg_masses, nothing)
        end
        next_free_state_index += dof(anext);
        push!(alast.children, anext); 
        alast = anext;
    end
    body = NPendulum(a0, -1);

    dof(body, force=true);

    return body;
end



# Attach one structure to another
function join_bodies!(parent_body::AbstractArticulatedBody, parent_articulation::Articulation, child_body::AbstractArticulatedBody)
    for child_articulation in child_body.articulation_zero.children
        join_bodies!(parent_body, parent_articulation, child_articulation)
    end
end
function join_bodies!(parent_body::AbstractArticulatedBody, parent_articulation::Articulation, child_root::Articulation)
    push!(parent_articulation.children, child_root)
    child_root.parent = parent_articulation;
    refresh_all_indices!(parent_body);
end

# Attach one structure to the tip of another
function attach_to_tip!(parent_body::AbstractArticulatedBody, child_body::AbstractArticulatedBody)
    tips = get_tips(parent_body);
    if length(tips) != 1
        throw(ExceptionError("attach_to_tip only implemented for single-tip structures"))
    end
    join_bodies!(parent_body, tips[1], child_body)
end

# Must be carried out after doing operations which alter a bodies' dof.
function refresh_all_indices!(body::AbstractArticulatedBody)
    forward_recurse(body, refresh_all_indices!, (1, 1));
    dof(body, force=true);
    return;
end
function refresh_all_indices!(a::Articulation, next_indices::Tuple{<: Integer, <: Integer})
    next_body_number, next_state_index = next_indices;
    a.body_number = next_body_number; next_body_number += 1;
    a.state_indices = next_state_index:(next_state_index + dof(a) - 1); next_state_index += dof(a);
    return (next_body_number, next_state_index);
end

# Discretization calculations for the lengths of segments in a filament.
function even_discretization(total_length::Number, i::Integer, n::Integer)
    return total_length / n;
end
function gauss_lobatto_discretization(total_length::Number, i::Integer, n::Integer)
    delta_s = 0.5((1-cos(pi*i/n)) - (1-cos(pi*(i-1)/n)));
    return total_length * delta_s;
end
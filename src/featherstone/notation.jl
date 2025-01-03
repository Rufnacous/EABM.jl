
dof(a::Articulation) = dof(a.joint_type);
function dof(b::AbstractArticulatedBody; force::Bool=false)
    if (b.dofs == -1) || (force)
        b.dofs = forward_recurse(b, (a::Articulation, count::Integer) -> count+dof(a), 0);
    end
    return b.dofs
end

function n_bodies(b::AbstractArticulatedBody)
    return forward_recurse(b, (a::Articulation, count::Integer) -> count+1, 0);
end

function Base.getindex(harness::StateHarness, i::Articulation)
    return harness.articulation_states[i.body_number + 1];
end
function Base.setindex!(harness::StateHarness, a_h::ArticulationHarness, i::Articulation)
    harness.articulation_states[i.body_number + 1] = a_h;
end
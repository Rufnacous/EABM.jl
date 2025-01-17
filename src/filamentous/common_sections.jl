
# Example second moment and articulation definitions for simple geometries.

function cylinder_second_moment_area(radius)
    Ix_Iy = pi * (radius ^ 4) / 4;
    return [Ix_Iy 0 0; 0 Ix_Iy 0 ; 0 0 0];
end
function rectangular_beam_second_moment_area(width, thickness)
    return [width*width*width*thickness/12 0 0; 0 thickness*thickness*thickness*width/12 0; 0 0 0];
end

RodArticulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Union{Articulation,Nothing}, length::Real, radius::Real, density::Real, stiffness::Real;
    curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], pregeometry::Vector{<:Real}=[0,0,0]
) = Articulation(
    number, next_free_state_index, joint, parent, length,
    cylinder_mass_and_inertia(density, length, radius)..., curvature, [0,0,length], pregeometry, [0,0,0.5length],
    FilamentProperties(LinearElasticProperties(stiffness * cylinder_second_moment_area(radius) / length))
);


mutable struct Rod <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function Rod(rod_length::Number, rod_radius::Number, rod_density::Number, rod_stiffness::Number, n::Integer, joints::JointType)
    a0 = articulation_zero();
    alast =  a0; anext = nothing;

    next_free_state_index = 1;

    for ni = 1:n
        anext = RodArticulation(ni, next_free_state_index, joints, alast, rod_length/n, rod_radius, rod_density, rod_stiffness)
        next_free_state_index += dof(anext);
        push!(alast.children, anext); anext.parent = alast;
        alast = anext;
    end
    body = Rod(a0, -1);

    dof(body, force=true);

    return body;
end
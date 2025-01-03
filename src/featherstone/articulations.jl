λ(a::Articulation) = a.parent;

articulation_zero() = Articulation(
    0, 0:-1,
    ArticulationZeroJoint(),
    nothing, Articulation[],
    0, 0, zeros(6,6),
    [0,0,0], [0,0,0], [0,0,0],
    Matrix{Float64}(I(6)), Matrix{Float64}(I(6)),
    nothing
);


Articulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Union{Articulation,Nothing}, length::Real, mass::Real, inertia_matrix_about_center::Matrix{<:Real},
    curvature::Matrix{<:Real},
    child_anchor::Vector{<:Real}, pregeometry::Vector{<:Real}, mass_center::Vector{<:Real},
    properties
) = Articulation(
    number, next_free_state_index:(next_free_state_index+dof(joint)-1),
    joint, parent, Articulation[],
    length, mass, inertia_matrix_about_origin(mass, inertia_matrix_about_center, mass_center),
    pregeometry, child_anchor, mass_center, diag_spatial(curvature)*xlt(child_anchor), xlt(pregeometry),
    properties
);

function get_articulation(body::AbstractArticulatedBody, by_number::Number)
    return forward_recurse(body, (candidate::Articulation, storage::Union{Nothing, Articulation}) -> ifelse(candidate.body_number == by_number, candidate, storage), nothing );
end
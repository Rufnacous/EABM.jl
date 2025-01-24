
# Forces and torques are callable.
function (force::ExternalForce)(a::Articulation, i::ArticulationHarness, t::Real)
    # This transform converts the force from global to local terms for the ABA,
    # and applies the force to the mass center.
    return xlt(a.mass_center)' * i.X0 * [0,0,0,force.force(a,i,t)...];
end
(torque::InternalTorque)(args...) = torque.torque(args...);

# Overloading (+) to allow force and torque addition.
(+)(a::AbstractExternalForce, b::AbstractExternalForce) = CompoundExternalForces(a,b);
function (forces::CompoundExternalForces)(a::Articulation, i::ArticulationHarness, t::Real)
    return forces.a(a,i,t) .+ forces.b(a,i,t)
end

(+)(a::AbstractInternalTorque, b::AbstractInternalTorque) = CompoundInternalTorques(a,b);
function (torques::CompoundInternalTorques)(a::Articulation, i::ArticulationHarness, t::Real)
    return torques.a(a,i,t) .+ torques.b(a,i,t)
end

# Some easy example forces
force_none() = ExternalForce((args...; kwargs...) -> [0,0,0]);
torque_none() = InternalTorque((a::Articulation, args...; kwargs...) -> zeros(dof(a)));

torque_restorative(;k::Number=1) = InternalTorque((a::Articulation, i::ArticulationHarness, t::Real) -> -k .* i.q);
torque_damping(;c::Number=1) = InternalTorque((a::Articulation, i::ArticulationHarness, t::Real) -> -c .* i.q_dt);


force_gravity(;g::Number=9.81, downwards::Vector{<:Number}=[0,0,-1]) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> a.mass*g .* downwards
);


# Calculate reaction force at root

function get_reaction_force(b::AbstractArticulatedBody, q::Vector{<:Number}, force::AbstractExternalForce, torque::AbstractInternalTorque, t::Number; alg::ArticulatedBodyAlgorithm=featherstones_algorithm)
    s = StateHarness(Float64, b);
    alg(b, s, q, t, force, torque);

    return inv(s[b.articulation_zero.children[1]].X0) * s[b.articulation_zero.children[1]].pA;
end
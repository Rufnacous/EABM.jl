
function frequencies(body::AbstractArticulatedBody; dynamics_algorithm::ArticulatedBodyAlgorithm=featherstones_algorithm)
    return frequencies(body, force_none(), torque_elastic(), dynamics_algorithm=dynamics_algorithm);
end
function frequencies(body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque; dynamics_algorithm::ArticulatedBodyAlgorithm=featherstones_algorithm)
    j = jacobian( make_linear_analysis_problem(dynamics_algorithm, body, force, torque) , zeros(dof(body)) );

    F = eigen(j);
    eigenmodes = real.(F.vectors);
    eigenfrequencies = sqrt.(abs.(real.(F.values)));

    return eigenfrequencies, eigenmodes;
end

function linear_analysis_problem(u::Vector{<:Real}, body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque, dynamics_algorithm::ArticulatedBodyAlgorithm)
    h = StateHarness(typeof(u).parameters[1], body)
    return dynamics_algorithm(body, h, [u..., zeros(dof(body))...], 0, force, torque);
end

function make_linear_analysis_problem(dynamics_algorithm::ArticulatedBodyAlgorithm, body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque)
    return let dynamics_algorithm = dynamics_algorithm, body = body, force = force, torque = torque
        (u) -> linear_analysis_problem(u, body, force, torque, dynamics_algorithm);
    end
end
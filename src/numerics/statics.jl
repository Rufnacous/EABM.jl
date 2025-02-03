# Basic statics procedure with newton raphson iteration.
# Find q.. (accelerations) = 0.

function static_posture(body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque; 
        solver=newton_raphson(0.01), repeat_condition=standard_repeat_condition)
    float_state_harness = StateHarness(Float64, body);
    dual_state_harness = StateHarness(Float64, body);

    J = (x,p) -> jacobian((u) -> statics_problem(u,p), x);
    
    initial = zeros(dof(body));
    return solver(statics_problem, J, initial, repeat_condition, (body, (float_state_harness, dual_harness), force, torque));
end




function newton_raphson(relaxation::Number)
    function solver(f, J, x::Vector{<: Real}, repeat_condition, p)
        i = 0;
        x_next = x .- (inv(J(x, p)) * f(x, p));
        while repeat_condition(i, x, x_next, f, J)
            x = x .+ relaxation * (x_next .- x);
            x_next = x .- (inv(J(x)) * f(x));
            i += 1;
        end
        return x;
    end
end

function standard_repeat_condition(i::Integer, x::Vector{<: Real}, x_next::Vector{<: Real}, f, J)
    return maximum(abs.(x_next .- x)) > 0.001
end





function statics_problem(u,p,t)
    ab_torque_calculator(p[1], pick_harness(typeof(u).parameters[1], p[2]...), u, 0, p[3], p[4])
end
function statics_problem(u,p)
    statics_problem(u,p,0)
end

function aba_pass2_statics!(
    a::Articulation, i::ArticulationHarness,
    aλ::Articulation, iλ::ArticulationHarness,
    t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

    i.U = i.IA * i.S;
    i.D = i.S' * i.U;
    i.D⁻¹ = inv(i.D);
    stiffnesses = diag(i.S[1:3,:]' * a.properties.elastics.spatial_rigidity * i.S[1:3,:]);
    i.u = (τ(a,i,t) - (i.S' * i.pA)) ./ stiffnesses
    # decreases problem stiffness + improves ease of use.
    # u is returned by get_torque but corresponds to a linear estimate of the
    # displaced angle.
    
    if aλ.body_number == 0
        return
    end
    Ia = i.IA
        # - (i.U * i.D⁻¹ * (i.U)');
    pa = i.pA
        # + (Ia * i.c) + (i.U * i.D⁻¹ * i.u);
    # removed to decrease stiffness, as unecessary to solution.

    iλ.IA += i.λX⁻ᵀ * Ia * i.Xλ;
    iλ.pA += i.λX⁻ᵀ * pa;
    return

end

# ab_torque_calculator is an ArticulatedBodyAlgorithm for use in statics solution.
ab_torque_calculator = ArticulatedBodyAlgorithm(
    [[ ((a,s,t,fx,τ) -> step(a, s[a], λ(a), s[λ(a)], t, fx, τ), action)
        for (step, action) in 
            [ (aba_pass1!, :forward),
            (aba_pass2_statics!, :backward)
            ] ]..., (get_torque, :return)]
    );

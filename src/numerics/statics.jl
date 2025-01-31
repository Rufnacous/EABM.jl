# Basic statics procedure with newton raphson iteration.
# Find q.. (accelerations) = 0.

function static_posture(body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque; 
        solver=newton_raphson(0.01), repeat_condition=standard_repeat_condition)
    float_state_harness = StateHarness(Float64, body);
    dual_state_harness = StateHarness(Float64, body);
    f = make_statics_problem(body, force, torque, float_state_harness);

    Jfx = make_statics_problem(body, force, torque, dual_state_harness);
    J(x) = jacobian(u -> Jfx(u), x);
    
    initial = zeros(dof(body));
    return solver(f, J, initial, repeat_condition);
end




function newton_raphson(relaxation::Number)
    function solver(f, J, x::Vector{<: Real}, repeat_condition)
        i = 0;
        x_next = x .- (inv(J(x)) * f(x));
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







function make_statics_problem(body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque, state_harness::StateHarness)
    return let body = body, force = force, torque = torque, state_harness = state_harness
        (du,u) -> ab_torque_calculator(body, state_harness, u, 0, force, torque,du)
    end
end

function statics_problem(du,u,p,t)
    ab_torque_calculator(p[1], p[2], u, 0, p[3], p[4], du)
end

function aba_pass2_statics!(
    a::Articulation, i::ArticulationHarness,
    aλ::Articulation, iλ::ArticulationHarness,
    t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

    mul!(i.U, i.IA, i.S);
    mul!(i.D, i.S', i.U);
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
            ] ]..., (get_torque!, :return)]
    );

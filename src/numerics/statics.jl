
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
        calculator = ab_torque_calculator(body);
        (u) -> calculator(body, state_harness, u, 0, force, torque)
    end
end

function aba_pass2_statics(stiffnesses)

    function aba_pass2_statics!(
        a::Articulation, i::ArticulationHarness,
        aλ::Articulation, iλ::ArticulationHarness,
        t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

        i.U = i.IA * i.S;
        i.D = i.S' * i.U;
        i.D⁻¹ = inv(i.D);
        i.u = (τ(a,i,t) - (i.S' * i.pA)) ./ stiffnesses[a.state_indices];
        
        if aλ.body_number == 0
            return
        end
        Ia = i.IA # - (i.U * i.D⁻¹ * (i.U)');
        pa = i.pA # + (Ia * i.c) + (i.U * i.D⁻¹ * i.u);

        iλ.IA += i.λX⁻ᵀ * Ia * i.Xλ;
        iλ.pA += i.λX⁻ᵀ * pa;
        return
    end

    return aba_pass2_statics!;
end

ab_torque_calculator(body::AbstractArticulatedBody) = ArticulatedBodyAlgorithm(
    [[ ((a,s,t,fx,τ) -> step(a, s[a], λ(a), s[λ(a)], t, fx, τ), action)
        for (step, action) in 
            [ (aba_pass1!, :forward),
            (aba_pass2_statics(get_stiffnesses(body)), :backward)
            ] ]..., (get_torque, :return)]
    );


# An ArticulatedBodyAlgorithm is a callable struct.

function (algorithm::ArticulatedBodyAlgorithm)(articulated_body::AbstractArticulatedBody,
    # this caller creates a StateHarness so you dont have to.
        generalized_state::Vector{<: Real}, time::Real,
        force::AbstractExternalForce, torque::AbstractInternalTorque)
    return algorithm(articulated_body, StateHarness(Float64, articulated_body), generalized_state, time, force, torque);
end
function (algorithm::ArticulatedBodyAlgorithm)(articulated_body::AbstractArticulatedBody,
    harness::StateHarness, generalized_state::Vector{<: Real}, time::Real,
    force::AbstractExternalForce, torque::AbstractInternalTorque)

    # Always set_state! before calling ABA steps.
    set_state!(articulated_body, harness, generalized_state);

    # Each step is either the final step which returns data, or is forward/backward recursed.
    for (step, action) in algorithm.passes
        if (action == :return)
            return step(articulated_body, harness);
        elseif (action == :forward)
            forward_recurse_iterative!(articulated_body, step, harness, time, force, torque);
        elseif (action == :backward)
            backward_recurse_iterative!(articulated_body, step, harness, time, force, torque);
        end
    end
    
end

# The Articulated Body Algortihm. See Page 132 of Rigid Body Dynamics Algorithms by Roy Featherstone (2008).

function aba_pass1!(
    a::Articulation, i::ArticulationHarness,
    aλ::Articulation, iλ::ArticulationHarness,
    t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )
    
    i.IA = a.I;
    i.pA = i.v ⨳ (a.I * i.v) - fx(a,i,t);
    return
end

function aba_pass2!(
    a::Articulation, i::ArticulationHarness,
    aλ::Articulation, iλ::ArticulationHarness,
    t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

    i.U = i.IA * i.S;
    i.D = i.S' * i.U;
    i.D⁻¹ = inv(i.D);
    i.u = τ(a,i,t) - (i.S' * i.pA);
    
    if aλ.body_number == 0
        return
    end
    Ia = i.IA - (i.U * i.D⁻¹ * (i.U)');
    pa = i.pA + (Ia * i.c) + (i.U * i.D⁻¹ * i.u);

    iλ.IA += i.λX⁻ᵀ * Ia * i.Xλ;
    iλ.pA += i.λX⁻ᵀ * pa;
    return
end

function aba_pass3!(
    a::Articulation, i::ArticulationHarness,
    aλ::Articulation, iλ::ArticulationHarness,
    t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )


    a_prime = (i.Xλ * iλ.a) + i.c;
    i.q_d²t = i.D⁻¹ * (i.u - (i.U' * a_prime));
    i.a = a_prime + (i.S * i.q_d²t);

end

# Composes these steps into an ArticulatedBodyAlgorithm.
featherstones_algorithm = ArticulatedBodyAlgorithm(
    [[ ((a,s,t,fx,τ) -> step(a, s[a], λ(a), s[λ(a)], t, fx, τ), action)
        for (step, action) in 
            [ (aba_pass1!, :forward),
            (aba_pass2!, :backward),
            (aba_pass3!, :forward),
            ] ]..., (get_acceleration, :return)]
    );



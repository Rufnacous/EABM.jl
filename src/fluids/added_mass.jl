

function featherstones_with_added_mass(fluid_density::Number, flow::Function)
    function added_mass_aba_pass1!(
        a::Articulation, i::ArticulationHarness,
        aλ::Articulation, iλ::ArticulationHarness,
        t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

        If, Sf = get_added_inertia_properties(a.properties.ffi, i, fluid_density);
        Uf = If * Sf; Df = Sf' * Uf; D⁻¹f = inv(Df);
        Iaf = If - (Uf * D⁻¹f * Uf');

        vJf = Sf*(Sf'* ([0,0,0,flow(i.p, t)...] .- i.v));
        vf = i.v + vJf;

        dudt = flow_acceleration(flow, i.p, t);
        pAf = Iaf * ( (vf ⨱ vJf) - (i.X0 * [0,0,0, (-dudt)...]) );

        i.IA = a.I #+ Iaf;
        i.pA = i.v ⨳ (a.I * i.v) - fx(a,i,t) + pAf;
        return
    end

    return ArticulatedBodyAlgorithm(
    [[ ((a,s,t,fx,τ) -> step(a, s[a], λ(a), s[λ(a)], t, fx, τ), action)
        for (step, action) in 
            [ (added_mass_aba_pass1!, :forward),
            (aba_pass2!, :backward),
            (aba_pass3!, :forward),
            ] ]..., (get_acceleration, :return)]
    );
end
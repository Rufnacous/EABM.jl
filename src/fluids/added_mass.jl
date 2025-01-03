

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

        i.IA = a.I + Iaf;
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


function force_further_added_mass_for_elongated_bodies(fluid_density::Number, flow_func::Function)
    function f_ma_elongated(a::Articulation, i::ArticulationHarness, t::Real)

        ma = a.properties.ffi.added_fluid_volume * fluid_density;
        
        u_global_absolute = [0,0,0, flow_func(i.p, t)...];
        u_local_absolute = inv(xlt(a.child_anchor)) * inv(i.X0') * u_global_absolute;
        u_local_relative = i.v - u_local_absolute;

        u_normal = u_local_relative[4];
        u_tangential = u_local_relative[6];
        theta_dot = i.V[2];
        curvature = 1 / (a.length / i.q[1]);
        
        # (-2theta_dot*u_tangential) + 
        force = - ma * ((curvature*((u_tangential^2)-0.5(u_normal^2)))) * i.X0[2,1:3];

        return ifelse(i.q[1] == 0, [0,0,0], force);
    end
    return ExternalForce(f_ma_elongated);
end
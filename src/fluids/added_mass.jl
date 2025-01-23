
# featherstones_with_added_mass is an ArticulatedBodyAlgorithm,
# as it makes modifications to the featherstones_algorithm step 1.
# This is described in the journal paper for the EABM. Section 2.3.2
function featherstones_with_added_mass(fluid_density::Number, flow::Function)
    function added_mass_aba_pass1!(
        a::Articulation, i::ArticulationHarness,
        aλ::Articulation, iλ::ArticulationHarness,
        t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

        If, Sf = get_added_inertia_properties(a.properties.ffi, i, fluid_density);
        Uf = If * Sf; Df = Sf' * Uf; D⁻¹f = inv(Df);
        Iaf = If - (Uf * D⁻¹f * Uf');
        # display(If)
        # # display(Sf)
        # display(Uf * D⁻¹f * Uf')
        # display(Iaf)
        # readline()

        U_global_absolute = [0,0,0,flow(i.p, t)...];
        U_local_absolute = i.X0 * U_global_absolute;
        U_local_relative = U_local_absolute .- i.v;

        vJf = Sf * Sf' * U_local_relative;
        vf = i.v + vJf;

        Du_global = flow_acceleration(flow, i.p, t);
        Du_local = i.X0 * [0,0,0, Du_global...];
        pAf = Iaf * ( (vf ⨱ vJf) + (vf ⨳ vf) - Du_local );

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

# See Leclercq and de Langre's 2018 paper "Reconfiguration of elastic blades in oscillatory flow"
function force_further_added_mass_for_elongated_bodies(fluid_density::Number, flow_func::Function)
    function f_ma_elongated(a::Articulation, i::ArticulationHarness, t::Real)

        ma = a.properties.ffi.added_fluid_volume * fluid_density;
        
        U_global_absolute = [0,0,0,flow_func(i.p, t)...];
        U_local_absolute = i.X0 * U_global_absolute;
        U_local_relative = U_local_absolute - i.v

        u_normal = - U_local_relative[4]; # of solid relative to fluid
        u_tangential = - U_local_relative[6]; # of solid relative to fluid
        
        curvature = 1 / (a.length / i.q[1]);
        force = - ma * ((curvature*((u_tangential^2)-0.5(u_normal^2)))) * i.X0[1,1:3];

        return ifelse(i.q[1] == 0, [0,0,0], force);
    end
    return ExternalForce(f_ma_elongated);
end


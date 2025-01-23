using EABM
using Test
import Plots.default, Plots.plot, Plots.plot!, Plots.cgrad, Plots.@layout, Plots.palette, Plots.savefig;
import Printf.@sprintf;

@testset "EABM.jl linear" begin
    
    mkpath("./output");

    n_cycles=2;
    n = 10; rod_length = 0.5; radius = 0.01; density = 10; stiffness = 1e3;
    
    theoretical_freqs = (([4.694, 1.875] .^ 2) ./ (rod_length^2)) .* sqrt(stiffness * EABM.cylinder_second_moment_area(radius)[1] / (density * (pi * radius^2))) 
    println("Theoretical frequencies = ", theoretical_freqs)

    T = n_cycles * 2pi / theoretical_freqs[end];

    body = Rod(rod_length, radius, density, stiffness, n, RotaryJoint(:x));
    force = force_none(); torque = torque_elastic();

    predicted_freqs, modes = frequencies(body, force, torque);
    println("Predicted frequencies =   ",predicted_freqs[end-1:end] );

    iq = zeros(body);
    iq[1:dof(body)] .= 0.1 * modes[:,end];
    sol = simulate(body, force, torque, T, initcond = iq, dt_modify = 0.1, integrator = :stable);

    p1 = plot(sol, idxs=1:n, xlabel="time / s", ylabel="qⁱ / radians", title="Oscillation Time History", label="");
    cols = cgrad(:roma);
    p2 = plot([],[],label="", xlims=(-0.25,0.25), ylims=(0,0.5), title="Oscillation Frames", xlabel="y axis /m", ylabel="z axis / m");
    ts = LinRange(((n_cycles-1)/n_cycles)*sol.t[end], sol.t[end], 50);
    for ti in eachindex(ts)
        t = ts[ti];
        pos = get_position(body, sol(t));
        plot!(pos[2,:], pos[3,:], label="", color=cols[ti  / length(ts)])
    end

    l = @layout [a b];
    plot(p1, p2, layout=l);
    
    savefig("./output/linear.png");

end


@testset "EABM.jl statics" begin
    
    fluid_density = 1000;
    torque = torque_elastic();

    materialA(L, n) = EABM.FluidborneStrip(L, 0.015, 0.006, 1190, 3200e6, n, RotaryJoint(:y), Cd=1.95);
    materialB(L, n) = EABM.FluidborneStrip(L, 0.015, 0.0013, 1245, 1620e6, n, RotaryJoint(:y), Cd=1.95);
    materialC(L, n) = EABM.FluidborneStrip(L, 0.015, 0.0004, 1055, 904e6, n, RotaryJoint(:y), Cd=1.95, curvature=EABM.rotate_y(L/(n * 0.436)));

    surrogate1 = materialA(0.3, 30);
    
    surrogate2 = materialA(0.1, 10);
    EABM.attach_to_tip!(surrogate2, materialB(0.1, 10));
    EABM.attach_to_tip!(surrogate2, materialC(0.1, 10));
    
    surrogate3 = materialA(0.1, 10);
    EABM.attach_to_tip!(surrogate3, materialC(0.2, 20));
    
    surrogate4 = materialC(0.3, 30);

    test_subjects = [ surrogate1, surrogate2, surrogate3, surrogate4 ];
    flows = [(0.0031, 0.0499), (0.0066, 0.1078), (0.0084, 0.1979), (0.0215, 0.3197)];

    plot([],[],label="", xlims=(-0.01, 1.2), ylims=(0,0.5),title="Heterogeneous Plastic Strips Under Varied Fluid Loads", aspect_ratio=:equal);

    println("Starting statics... ")

    for body_i in eachindex(test_subjects)
        body = test_subjects[body_i];
        println(" Solving for body ", body_i)
        for flow_i in eachindex(flows)
            println("             flow ", flow_i)
            flow = flows[flow_i];
            
            flow_func = (x,t)->[flow[1]*log(max(0.0001,x[3])) + flow[2],0,0];
            
            force = EABM.force_buoyancy(fluid_density) + 
                EABM.force_gravity() +
                EABM.force_drag(fluid_density, flow_func)

            static_result = static_posture(body, force, torque);


            pos = get_position(body, static_result);
            plot!(pos[1,:] .+ 0.3*(flow_i-1), pos[3,:],
                label=ifelse(flow_i==1,@sprintf("Surrogate %d",body_i),""),
                color=palette(:tab10)[body_i], linewidth=2)
        end
    end

    mkpath("./output");
    savefig("./output/statics.png");
    
end




@testset "EABM.jl dynamics" begin

    function wave_velocity(pos, t)
        k = 1.9287; ramptime = 10;
        z = pos[3];
        x = pos[1];

        U_A_M = [ 0.0694  0.1175  0.0762  0.0414];
        U_A_C = [ 0.1582  0.0462  0.0075  0.0001];
        U_P   = [-1.5072 -2.3098 -2.7410  1.2825];
        W_A_M = [ 0.2950  0.1520  0.0887  0.0116];
        W_A_C = [ 0.0129  0.0118  0.0032  0.0014];
        W_P   = [ 0.1719 -0.5177 -1.3489 -1.2099];

        U_A = U_A_M*z + U_A_C;
        W_A = W_A_M*z + W_A_C;

        U = 0; W = 0;
        for harmonic = 1:4
            U += U_A[harmonic] * cos(pi*harmonic*(t) + U_P[harmonic] - harmonic*k*x);
            W += W_A[harmonic] * cos(pi*harmonic*(t) + W_P[harmonic] - harmonic*k*x);
        end

        if t < 0
            U = 0;
            W = 0;
        elseif t < ramptime
            U *= (1-cos(t*pi/ramptime))/2;
            W *= (1-cos(t*pi/ramptime))/2;
        end
        wave_vel = [U, 0, W];

        return wave_vel;
    end

    n = 50; T=30;

    rod_length = 0.2; width = 0.02; thickness = 0.0019; density = 670; stiffness = 0.5e6;
    
    # See Zhu et al. 2020 "Mechanisms for the Asymmetric Motion..." for Cm calculation
    KC = 0.206*2/width;
    Cm = min(ifelse(KC < 20, 1+0.35(KC^2/3), 1+0.15(KC^2/3)), 1+((KC-18)^2)/49);

    body = FluidborneStrip(rod_length, width, thickness, density, stiffness, n, RotaryJoint(:y), Cm = Cm,
        discretization = EABM.gauss_lobatto_discretization);

    fluid_density = 1000;
    flow_func = wave_velocity;

    force = force_gravity() + force_buoyancy(fluid_density) + 
        force_drag(fluid_density, flow_func) +
        force_skin_friction(fluid_density, flow_func) + 
        force_virtual_buoyancy(fluid_density, flow_func) + 
        force_further_added_mass_for_elongated_bodies(fluid_density, flow_func)

    torque = torque_elastic();
    
    predicted_freqs, modes = frequencies(body, force_none(), torque,
        dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func)
    );
    println("Predicted frequencies =   ",predicted_freqs[end-1:end] );
    println("                           [4.69, 0.36]")

    iq = zeros(body);


    sol = simulate(body, force, torque, T, initcond = zeros(body), 
        dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func),
        integrator=:approx, checkpoints = collect(LinRange(0,1,21))
        );


    cols = cgrad(:roma);
    plot([],[],label="", xlims=(-0.05,0.15), ylims=(0,0.2),aspect_ratio=:equal, xticks=[-0.05, 0, 0.05, 0.1, 0.15], yticks=[0,0.05,0.1,0.15,0.2],title="Silicon Strip in Wave Tank");
    ts = LinRange(sol.t[end] - 2, sol.t[end], 50);
    for ti in eachindex(ts)
        t = ts[ti];
        pos = get_position(body, sol(t));
        plot!(pos[1,:], pos[3,:], label="", color=cols[ti  / length(ts)])
    end
    mkpath("./output");
    savefig("./output/dynamics.png");




end

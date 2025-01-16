using EABM
using Test
import Plots.default, Plots.plot, Plots.plot!, Plots.cgrad, Plots.@layout, Plots.palette;
import Printf.@sprintf;

@testset "EABM.jl linear" begin

    n_cycles=2;
    n = 10; rod_length = 0.5; radius = 0.01; density = 10; stiffness = 1e3;
    
    theoretical_freqs = (([4.694, 1.875] .^ 2) ./ (rod_length^2)) .* sqrt(stiffness * EABM.cylinder_second_moment_area(radius) / (density * (pi * radius^2))) 
    println("Thheoretical frequencies = ", theoretical_freqs)

    T = n_cycles * 2pi / theoretical_freqs[end];

    body = Rod(rod_length, radius, density, stiffness, n, RotaryJoint(:x));
    force = force_none(); torque = torque_elastic();

    predicted_freqs, modes = frequencies(body, force, torque);
    println("Predicted frequencies =   ",predicted_freqs[end-1:end] );

    iq = zeros(body);
    iq[1:dof(body)] .= 0.1 * modes[:,end];
    sol = simulate(body, force, torque, T, initcond = iq, dt_modify = 0.1, integrator = :stable);

    default(show=true);

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
    println(" See results ")

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

    plot([],[],label="", xlims=(-0.01, 1.2), aspect_ratio=:equal);

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

            default(show=true);

            pos = get_position(body, static_result);
            plot!(pos[1,:] .+ 0.3*(flow_i-1), pos[3,:],
                label=ifelse(flow_i==1,@sprintf("Body %d",body_i),""),
                color=palette(:tab10)[body_i], linewidth=2)
        end
    end
    
end

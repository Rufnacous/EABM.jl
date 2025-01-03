using EABM
using Test
import Plots.default, Plots.plot, Plots.plot!, Plots.cgrad, Plots.@layout;

@testset "EABM.jl linear" begin

    n_cycles=2;
    n = 10; rod_length = 0.5; radius = 0.01; density = 10; stiffness = 1e3;
    
    theoretical_freqs = (([4.694, 1.875] .^ 2) ./ (rod_length^2)) .* sqrt(stiffness * EABM.cylinder_second_moment_area(radius) / (density * (pi * radius^2))) 
    println("Pheoretical frequencies = ", theoretical_freqs)

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
    sleep(10)

end

module EABM

    using LinearAlgebra;
    import Base.+;
    import Printf.@sprintf;
    import DifferentialEquations.ODEProblem, DifferentialEquations.solve, DifferentialEquations.SSPRK22, DifferentialEquations.TRBDF2;
    import ForwardDiff.Dual, ForwardDiff.jacobian;
    import Plots.default, Plots.plot, Plots.plot!, Plots.scatter, Plots.scatter!;

    include("featherstone/types.jl");
    include("featherstone/common_geometry.jl");
    include("featherstone/common_bodies.jl");
    include("featherstone/common_state.jl");
    include("featherstone/recursion.jl");
    include("featherstone/notation.jl");
    include("featherstone/spatial_algebra.jl");
    include("featherstone/articulations.jl");
    include("featherstone/joints.jl");
    include("featherstone/algorithm.jl");
    include("featherstone/driving.jl");

    include("numerics/linear.jl");
    include("numerics/statics.jl");
    include("numerics/simulation.jl");
    
    include("filamentous/types.jl");
    include("filamentous/common_sections.jl");
    include("filamentous/elastics.jl");

    include("fluids/types.jl");
    include("fluids/common_vectors.jl");
    include("fluids/common_subjects.jl");
    include("fluids/basic_forces.jl");
    include("fluids/added_mass.jl");

    export frequencies, static_solution, simulate;
    export NPendulum, Rod;

    function test(;T=10, dt_modify=0.5, n = 2, flowvel =0.1)
        length = 0.5; radius = 0.01; density = 10; stiffness = 1e6;

        body = FluidborneRod(length, radius, density, stiffness, n, RotaryJoint(:x));

        fluid_density = 1000;
        flow_func = (xyz, t) -> ([0,flowvel, 0] .* sin(50pi * t)) * (1 - cos(pi * max(0, min(t, 0.5T)) / (0.5T)));
        force = force_drag(fluid_density, flow_func);

        # force_gravity(downwards=[0,1,0]); # force_none();

        torque = torque_elastic(); #+ torque_damping(c=0.01);

        freqs, modes = frequencies(body, force, torque,
            dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func)
            );
        println(freqs);
        # display(modes);

        iq = zeros(body);
        # iq[1] = 0.1;

        sol = simulate(body, force, torque, T,
            initcond = iq, dt_modify = dt_modify,
            dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func)
            );


        default(show=true);
        plot(sol, idxs=1:n);
        print(" plotted sol. "); readline();

        plot([],[],label="", xlims=(-0.25,0.25), ylims=(0,0.5));
        for t in LinRange(0.5T, T, 20)[1:end-1]
            pos = get_position(body, sol(t));
            
            plot!(pos[2,:], pos[3,:], label="", color=:black)
        end


        return;

    end

    

    function test_branched(;T=10, dt_modify=0.5, n = 2, flowvel =0.1)
        limb_length = 0.25; radius = 0.01; density = 10; stiffness = 1e6;

        body = FluidborneRod(limb_length, radius, density, stiffness, n, RotaryJoint(:x));
        tip = get_articulation(body, n);

        branch_angle_1 = rotate_y(pi / 4);
        branch1 = FluidborneRod(limb_length, radius, density, stiffness, n, RotaryJoint(:x));
        branch1.articulation_zero.children[1].geometry_transform = diag_spatial(branch_angle_1)*xlt(branch1.articulation_zero.children[1].child_anchor);
        join_bodies!(body, tip, branch1)

        branch_angle_2 = rotate_y(-pi / 4);
        branch2 = FluidborneRod(limb_length, radius, density, stiffness, n, RotaryJoint(:x));
        branch2.articulation_zero.children[1].geometry_transform = diag_spatial(branch_angle_2)*xlt(branch2.articulation_zero.children[1].child_anchor);
        join_bodies!(body, tip, branch2)





        fluid_density = 1000;
        flow_func = (xyz, t) -> ([0,flowvel, 0] .* sin(50pi * t)) * (1 - cos(pi * max(0, min(t, 0.5T)) / (0.5T)));
        force = force_drag(fluid_density, flow_func);
        torque = torque_elastic();

        freqs, modes = frequencies(body, force, torque,
            # dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func)
            );
        println(freqs);

        iq = zeros(body);

        sol = simulate(body, force, torque, T,
            initcond = iq, dt_modify = dt_modify,
            # dynamics_algorithm = featherstones_with_added_mass(fluid_density, flow_func)
            );


        default(show=true);
        plot(sol, idxs=1:3n);
        print(" plotted sol. "); readline();

        # plot([],[],[],label="") #, xlims=(-0.25,0.25), ylims=(0,0.5));
        for t in LinRange(0.5T, T, 100)[1:end-1]
            pos = get_position(body, sol(t), with_retrace=true);
            plot( pos[2,:], pos[3,:], label="", color=:black, #pos[1,:],
            # xlims=(-0.25, 0.25), ylims=(-0.25, 0.25), zlims=(0, 0.5)
            );
            readline();
        end


        return;

    end

    
#     function test_speed()
#         length = 0.5; radius = 0.01; density = 1000; stiffness = 1e6;
        
#         force = force_none;
#         torque = torque_elastic;

#         ns = Int.(floor.(10 .^ LinRange(0, 3, 100)));
#         times = zeros(size(ns));
#         n_repeats = 100;
            
#         for ni in eachindex(ns)

#             n = ns[ni];

#             a1 = RodArticulation(1, 1, RotaryJoint(:x), nothing, length, radius, density, stiffness);
#             body = ArticulatedBody(a1);

#             for nii in 2:n
#                 a2 = RodArticulation(nii, nii, RotaryJoint(:x),   a1,    length, radius, density, stiffness);
#                 push!(a1.children, a2);

#                 a1 = a2;
#             end

#             println(dof(body, force=true));
#             harness = StateHarness(Float64, body); state = zeros(body);

#             for repeat = 1:n_repeats
#                 start = time();
#                 featherstones_algorithm(body, harness, state, 0, force, torque);
#                 times[ni] += time() - start;
#             end
#             times[ni] /= n_repeats;
            
#         end

#         plot(log.(10,ns), log.(10,times))

#         return;

#     end
end

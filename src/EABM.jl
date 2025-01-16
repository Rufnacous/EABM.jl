module EABM

    using LinearAlgebra;
    import Base.+;
    import Printf.@sprintf;
    import DifferentialEquations.ODEProblem, DifferentialEquations.solve, DifferentialEquations.SSPRK22, DifferentialEquations.TRBDF2, DifferentialEquations.DiscreteCallback;
    import ForwardDiff.Dual, ForwardDiff.jacobian;
    import Plots.default, Plots.plot, Plots.plot!, Plots.scatter, Plots.scatter!, Plots.cgrad, Plots.savefig;
    using CSV, Tables;

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
    export Rod, FluidborneRod, RotaryJoint;
    export force_none, torque_elastic, force_drag;
    export dof, n_bodies, get_positio;

    function test(;dt_modify=0.5, n = 2, T=8)
        rod_length = 0.2; width = 0.02; thickness = 0.0019; density = 670; stiffness = 0.5e6;

        body = FluidborneStrip(rod_length, width, thickness, density, stiffness, n, RotaryJoint(:x),
            discretization = gauss_lobatto_discretization);

        fluid_density = 1000;
        flow_func = asymmetric_wave_profile(0.039, 2pi / 2, 1.93, 0.3, 0.5T);
        force = force_drag(fluid_density, flow_func) + force_skin_friction(fluid_density, flow_func) + 
            force_gravity() + force_buoyancy(fluid_density) +
            force_further_added_mass_for_elongated_bodies(fluid_density, flow_func);
        torque = torque_elastic();
        
        predicted_freqs, modes = frequencies(body, force_none(), torque,
            dynamics_algorithm = featherstones_with_added_mass_and_virtual_buoyancy(fluid_density, flow_func)
        );
        println("Predicted frequencies =   ",predicted_freqs[end-1:end] );

        iq = zeros(body);

        sol = simulate(body, force, torque, T,
            initcond = iq, dt_modify = dt_modify,
            dynamics_algorithm = featherstones_with_added_mass_and_virtual_buoyancy(fluid_density, flow_func),
            integrator=:stable, checkpoints = collect(LinRange(0,1,101))
            );


        default(show=true);
        plot(sol, idxs=1:n);
        print(" plotted sol. "); readline();

        cols = cgrad(:roma);

        plot([],[],label="", xlims=(-0.125,0.125), ylims=(0,0.25));
        ts = LinRange(sol.t[end] - 2, sol.t[end], 50); #[2:end];
        for ti in eachindex(ts)
            t = ts[ti];
            pos = get_position(body, sol(t));
            plot!(pos[2,:], pos[3,:], label="", color=cols[ti  / length(ts)])
            readline()
        end



        return;

    end

    function SpeedTestRod(rod_length::Number, rod_radius::Number, rod_density::Number, rod_stiffness::Number, breadth::Integer, depth::Integer, joints::JointType)
        a0 = articulation_zero();

        next_free_state_index = 1;
        next_body = 1;
        a1 = RodArticulation(next_body, next_free_state_index, joints, a0, rod_length/depth, rod_radius, rod_density, rod_stiffness);
        next_free_state_index += dof(a1); next_body += 1;
        push!(a0.children, a1); a1.parent = a0;
    
        
        for bi = 1:breadth
            alast = a1;
            for di = 1:depth
                anext = RodArticulation(next_body, next_free_state_index, joints, alast, rod_length/depth, rod_radius, rod_density, rod_stiffness)
                next_free_state_index += dof(anext); next_body += 1;

                push!(alast.children, anext); anext.parent = alast;
                alast = anext;
            end
        end

        body = Rod(a0, -1);
        dof(body, force=true);
    
        return body;
    end

    
    function test_speed(;samples=10,repeatfor = 10.0, repeatminimum=20, ds = Int64.(floor.(10 .^ LinRange(0,2, 5))), bs = Int64.(floor.(10 .^ LinRange(0,2, 6))))
        rod_length = 0.2; radius = 0.0019; density = 670; stiffness = 0.5e6;
        
        force = force_none();
        torque = torque_elastic() + torque_damping() + InternalTorque((a::Articulation, i::ArticulationHarness, t::Real) -> ifelse(a.body_number==1,[3.4*sin(1.2*t)],[0]));
        alg = featherstones_algorithm;

        plot([],[],label="");

        dataout_mean = zeros(length(ds)+1, length(bs)+1);
        dataout_stdev = zeros(length(ds)+1, length(bs)+1);
        datasamples = zeros(length(ds), length(bs), samples);

        for bi in eachindex(bs)
            print(bi, " ")

            ts = zeros(length(ds));

            for di in eachindex(ds)
                print(di)
                body = SpeedTestRod(rod_length, radius, density, stiffness, bs[bi], ds[di], RotaryJoint(:x));

                s = StateHarness(Float64, body); u = zeros(body);

                for si in 1:samples
                    t = time()
                    for repeat in 1:repeatminimum
                    # n_repeats = 0; repeatedfor = 0;
                    # while (n_repeats < repeatminimum) #(repeatedfor < repeatfor) ||
                        alg(body, s, u, 0, force, torque);
    #                    repeatedfor = time() - t; n_repeats += 1;
                    end
                    repeatedfor = time() - t;
                    ts[di] = repeatedfor / repeatminimum;
                    # dataout[1+di, 1+bi] = ts[di];
                    datasamples[di,bi, si] = ts[di];
                end

                dataout_mean[1+di, 1+bi] = sum(datasamples[di,bi,:]) / samples;

                dataout_stdev[1+di, 1+bi] = sqrt(sum((datasamples[di,bi,:] .- dataout_mean[1+di, 1+bi]) .^ 2) / samples);

                print(@sprintf(" (%.2f,%.2f) ",1000*dataout_mean[1+di, 1+bi], 1000*dataout_stdev[1+di, 1+bi]))
                
                # println(" ",ts[di])

                dataout_mean[1+di, 1] = ds[di];
                dataout_stdev[1+di, 1] = ds[di];

            end
            dataout_mean[1, 1+bi] = bs[bi];
            dataout_stdev[1, 1+bi] = bs[bi];

            ns = ds .* bs[bi];

            q =  1 ./ (ts ./ ns);
            scatter!(log.(10, ns), log.(10, q), xlabel="log10 n", ylabel="log10 aba speed hz", label=@sprintf("b%d", bs[bi])) #, ylims=(1500,3500))
            println()

        end

        

        
        f = open("speedtest_mean.csv","w+");
        CSV.write(f, Tables.table(dataout_mean));
        close(f)
        f = open("speedtest_stdv.csv","w+");
        CSV.write(f, Tables.table(dataout_stdev));
        close(f)

        savefig("khz.png")
        return;

    end
        
    function asymmetric_wave_profile(aw, 𝜔, k, h, ramptime)
        function wave_velocity(pos, t)
            z = pos[3];
            x = pos[2];

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
            wave_vel = [0, U, W];

            return wave_vel;
        end
        return wave_velocity
    end


    function linear_wave_profile(aw, 𝜔, k, h, ramptime)
        function wave_velocity(pos, t)
            z = pos[3];
            x = pos[1];
    
            U = aw * 𝜔 * (cosh(k*z)/sinh(k*h)) * cos.(k*x .- 𝜔*t);
            W = aw * 𝜔 * (sinh(k*z)/sinh(k*h)) * sin.(k*x .- 𝜔*t);
    
            if t < ramptime
                wave_vel = [0, U, W] .* (0.5-0.5cos(t*pi/ramptime));
            else
                wave_vel = [0, U, W] ;
            end
    
            return wave_vel;
        end
        return wave_velocity
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

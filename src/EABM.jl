module EABM

    using LinearAlgebra;
    import Base.+;
    import Printf.@sprintf;
    import DifferentialEquations.ODEProblem, DifferentialEquations.solve, DifferentialEquations.SSPRK22, DifferentialEquations.TRBDF2, DifferentialEquations.DiscreteCallback;
    import ForwardDiff.Dual, ForwardDiff.jacobian;

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

    export frequencies, static_posture, simulate, featherstones_with_added_mass;
    export ArticulatedBody, Rod, FluidborneStrip, RotaryJoint, BendingJoint;
    export torque_elastic;
    export force_none, force_gravity, force_buoyancy, force_drag,
        force_skin_friction, force_virtual_buoyancy,
        force_further_added_mass_for_elongated_bodies;
    export dof, n_bodies, get_position;
    

    function memory_test(;n=1)

        body = NPendulum(n, 0.1, 0.1, RotaryJoint(:y), radius=0.01);
        harness = StateHarness(Float64, body);
        generalized_state = zeros(body);

        memallocd = @allocated set_state!(body, harness, generalized_state);
        println("Allocated by set_state! : ", memallocd);

        f = force_none()
        t = torque_none()
        al = body.articulation_zero;
        a = al.children[1];
        i = harness[a]; il = harness[al];
        
        memallocd = @allocated aba_pass1!( a,i,al,il, 0.0, f,t);
        println("Allocated by aba_pass1! : ", memallocd);
        
        memallocd = @allocated aba_pass2!( a,i,al,il, 0.0, f,t);
        println("Allocated by aba_pass2! : ", memallocd);
        
        memallocd = @allocated aba_pass3!( a,i,al,il, 0.0, f,t);
        println("Allocated by aba_pass3! : ", memallocd);

        
        memallocd = @allocated featherstones_algorithm(body, harness, generalized_state, 0.0, f, t)
        println("Allocated by featherstones_algorithm : ", memallocd);

        
        

        # aba_pass1!(
        #     a::Articulation, i::ArticulationHarness,
        #     aλ::Articulation, iλ::ArticulationHarness,
        #     t::Real, fx::AbstractExternalForce, τ::AbstractInternalTorque    )

        


        # a = body.articulation_zero.children[1];
        # memallocd = @allocated j_calc(a.joint_type, harness[a]);
        # println("Allocated by j_calc : ", memallocd);

        

        # E = zeros(3,3); p = zeros(3); q = [0]; qdt = [0];
        # memallocd = @allocated get_transformations!(a.joint_type, q,qdt, S,Ṡ,E,p);
        # println("Allocated by get_transformations! : ", memallocd);

        return

    end


end


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
    export Rod, FluidborneStrip, RotaryJoint, BendingJoint;
    export torque_elastic;
    export force_none, force_gravity, force_buoyancy, force_drag,
        force_skin_friction, force_virtual_buoyancy,
        force_further_added_mass_for_elongated_bodies;
    export dof, n_bodies, get_position;
    

end


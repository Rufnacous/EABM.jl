
# Basic setup of a vector cauchy problem wrapping the ABA.
# Integration operated by DifferentialEquations.jl.

# The state harness needs to be declared + allocated externally to the vector_cauchy_problem.
# If it is not, it will require allocating for every time vector_cauchy_problem is run, which
# should be avoided. Reallocating it outside the definition of vector_cauchy_problem means
# that two must be defined, one to be used for direct calculation, another for autodifferentiation.
# Multiple dispatch of vector_cauchy_problem then allows us to use the appropriate harness.
# When introduced, this step represented a 20% speed-up in integrations.

function simulate(body::AbstractArticulatedBody, force::AbstractExternalForce, torque::AbstractInternalTorque, simtime::Real;
    initcond::Vector{<:Real}=zeros(body), dt_modify::Real=1.0, dynamics_algorithm::ArticulatedBodyAlgorithm=featherstones_algorithm,
    integrator::Symbol=:approx, checkpoints::Vector{<:Real}=zeros(0))

    float_state_harness = StateHarness(Float64, body);
    dual_state_harness = StateHarness(Float64, body);
    
    checkpoints = checkpoints .* simtime;
    
    prob = ODEProblem(
        make_vector_cauchy_problem(dynamics_algorithm), initcond, (0, simtime), (body, float_state_harness, dual_state_harness, force, torque),
        tstops=checkpoints);

    freqs, shapes = frequencies(body, force, torque, dynamics_algorithm=dynamics_algorithm);
    dt = 0.5dt_modify / maximum(freqs);
    println("dt ",dt)

    starttime = time();
    condition(u, t, integrator) = any(t .== checkpoints);
    affect!(integrator) = println("T = ",@sprintf("%.2f : %.1f%% @ %.2f mins",integrator.t, 100integrator.t/simtime, (time() - starttime)/60));
    cb_checkpoints = DiscreteCallback(condition, affect!);

    if integrator == :approx
        return solve(prob, TRBDF2(autodiff=true), callback=cb_checkpoints);
    elseif integrator == :stable
        return solve(prob, SSPRK22(), dt = dt, callback=cb_checkpoints);
    else
        throw(ExceptionError("Unknown integrator type requested of EABM.jl."))
    end
    
end

function make_vector_cauchy_problem(dynamics_algorithm::ArticulatedBodyAlgorithm)
    return let dynamics_algorithm = dynamics_algorithm
        (args...) -> vector_cauchy_problem(args..., dynamics_algorithm);
    end
end

function vector_cauchy_problem(du::Vector{<:Number}, u::Vector{<:Number},
    p::Tuple{AbstractArticulatedBody, StateHarness, StateHarness, AbstractExternalForce, AbstractInternalTorque},
    t::Real, dynamics_algorithm::ArticulatedBodyAlgorithm)

    du[1:dof(p[1])] = u[1+dof(p[1]):end];
    du[1+dof(p[1]):end] = dynamics_algorithm(p[1], pick_harness(typeof(u).parameters[1], p[2], p[3]), u, t, p[4], p[5]);
end

function pick_harness(u_type::DataType, float_harness::StateHarness, dual_harness::StateHarness)
    if u_type <: Dual
        return dual_harness;
    end
    return float_harness;
end
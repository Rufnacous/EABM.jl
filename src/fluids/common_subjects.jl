
function frontal_area(ffi::StandardFFIProperties, i::ArticulationHarness, relative_u::Vector{<: Number})
    return ffi.frontal_area;
end
function skin_area(ffi::StandardFFIProperties, i::ArticulationHarness, relative_u::Vector{<: Number})
    return ffi.skin_area;
end
function get_added_inertia_properties(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number)
    vf = ffi.added_fluid_volume;
    mf = fluid_density * vf;
    return mf * ffi.added_inertia, ffi.fluid_freedom_subspace;
end
function get_virtual_buoyancy_inertia(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number)
    return fluid_density * ffi.volume * ffi.virtual_buoyancy;
end

###  FLUID-BORNE ROD

function FluidborneRodArticulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Union{Articulation,Nothing}, length::Real, radius::Real, density::Real, stiffness::Real;
    curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], pregeometry::Vector{<:Real}=[0,0,0])

    CM = 1;
    
    Sf = [ 
        0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0;
        1 0 0 0
    ];

    If_xyz = (2/5) * radius^2;
    added_inertia = inertia_matrix_about_origin(1, diagm([If_xyz, If_xyz, If_xyz]), [0,0,0]);
    virtual_buoyancy = inertia_matrix_about_origin(1, zeros(3,3), [0,0,0]);

    return Articulation(
        number, next_free_state_index, joint, parent, length,
        cylinder_mass_and_inertia(density, length, radius)..., curvature, [0,0,length], pregeometry, [0,0,0.5length],
        FluidborneFilamentProperties(
            AxisymmetricElasticProperties(stiffness * cylinder_second_moment_area(radius) / length),
            StandardFFIProperties(2.0, 0.01, length*2radius, length*2pi*radius, length*pi*radius^2, CM*(length * pi * (2radius)^2), added_inertia, Sf, virtual_buoyancy)
            )
    );
end

mutable struct FluidborneRod <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function FluidborneRod(rod_length::Number, rod_radius::Number, rod_density::Number, rod_stiffness::Number, n::Integer, joints::JointType)
    a0 = articulation_zero();
    alast =  a0; anext = nothing;

    next_free_state_index = 1;

    for ni = 1:n
        anext = FluidborneRodArticulation(ni, next_free_state_index, joints, alast, rod_length/n, rod_radius, rod_density, rod_stiffness)
        next_free_state_index += dof(anext);
        push!(alast.children, anext); anext.parent = alast;
        alast = anext;
    end
    body = Rod(a0, -1);

    dof(body, force=true);

    return body;
end

###  FLUID-BORNE STRIP

function FluidborneStripArticulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Union{Articulation,Nothing}, length::Real, width::Real, thickness::Real, density::Real, stiffness::Real;
    curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], pregeometry::Vector{<:Real}=[0,0,0])

    CM = 1;
    
    Sf = [ 
        0 1 0 0;
        0 0 1 0;
        0 0 0 1;
        0 0 0 0;
        0 0 0 0;
        1 0 0 0
    ];

    If_xyz = (2/5) * (width/2)^2;
    added_inertia = inertia_matrix_about_origin(1, diagm([If_xyz, If_xyz, If_xyz]), [0,0,0]);
    virtual_buoyancy = inertia_matrix_about_origin(1, zeros(3,3), [0,0,0]);

    return Articulation(
        number, next_free_state_index, joint, parent, length,
        rectangular_beam_mass_and_inertia(density, length, width, thickness)..., curvature, [0,0,length], pregeometry, [0,0,0.5length],
        FluidborneFilamentProperties(
            BendingElasticProperties((stiffness * rectangular_beam_second_moment_area(width, thickness)[[1,2]] ./ length)...),
            StandardFFIProperties(3.6, 0.014, length*width, length*(2width + 2thickness), length*width*thickness, CM*(length * pi * (width/2)^2), added_inertia, Sf, virtual_buoyancy)
            )
    );
end

mutable struct FluidborneStrip <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function FluidborneStrip(strip_length::Number, strip_width::Number, strip_thickness::Number, strip_density::Number, strip_stiffness::Number, n::Integer, joints::JointType; discretization::Function=even_discretization)
    a0 = articulation_zero();
    alast =  a0; anext = nothing;

    next_free_state_index = 1;

    for ni = 1:n
        segment_length = discretization(strip_length, ni, n);
        anext = FluidborneStripArticulation(ni, next_free_state_index, joints, alast, segment_length, strip_width, strip_thickness, strip_density, strip_stiffness)
        next_free_state_index += dof(anext);
        push!(alast.children, anext); anext.parent = alast;
        alast = anext;
    end
    body = Rod(a0, -1);

    dof(body, force=true);

    return body;
end

# Various functions that must be implemented for an FFIProperties type.

function frontal_area(ffi::StandardFFIProperties, i::ArticulationHarness, relative_u::Vector{<: Number})
    return ffi.frontal_area;
end
function skin_area(ffi::StandardFFIProperties, i::ArticulationHarness, relative_u::Vector{<: Number})
    return ffi.skin_area;
end
function get_added_inertia_properties(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number)
    mf = fluid_density * ffi.added_fluid_volume;
    return mf * ffi.added_inertia, ffi.fluid_freedom_subspace;
end
function get_virtual_buoyancy_inertia(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number)
    return fluid_density * ffi.volume * ffi.virtual_buoyancy;
end



###  FLUID-BORNE STRIP

function FluidborneStripArticulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Union{Articulation,Nothing}, length::Real, width::Real, thickness::Real, density::Real, stiffness::Real,
    Cd::Real, Cf::Real, Cm::Real;
    curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], pregeometry::Vector{<:Real}=[0,0,0])

    Sf = [ 
        0 0;
        0 0;
        0 1;
        0 0;
        0 0;
        1 0
    ];
    added_inertia = inertia_matrix_about_origin(1, diagm([1,1,1]) .* 0.4((length/2)^2), [0,0,0]);

    return Articulation(
        number, next_free_state_index, joint, parent, length,
        rectangular_beam_mass_and_inertia(density, length, width, thickness)...,
        curvature, [0,0,length], pregeometry, 
        [0,0,0.5length],
        FluidborneFilamentProperties(
            LinearElasticProperties((stiffness * rectangular_beam_second_moment_area(width, thickness) ./ length)),
            StandardFFIProperties(
                Cd, Cf,
                length*width, length*(2width + 2thickness),
                length*width*thickness,
                Cm*(length * pi * (width/2)^2),
                added_inertia, Sf)
            )
    );
end

mutable struct FluidborneStrip <: AbstractArticulatedBody
    articulation_zero::Articulation
    dofs::Integer
end
function FluidborneStrip(strip_length::Number, strip_width::Number,
        strip_thickness::Number, strip_density::Number, strip_stiffness::Number,
        n::Integer, joints::JointType; discretization::Function=even_discretization,
        curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], Cd::Real=3.6479, Cf::Real=0.014, Cm::Real=1)
    a0 = articulation_zero();
    alast =  a0; anext = nothing;

    next_free_state_index = 1;

    for ni = 1:n
        segment_length = discretization(strip_length, ni, n);
        anext = FluidborneStripArticulation(ni, next_free_state_index, joints, alast, segment_length, strip_width, strip_thickness, strip_density, strip_stiffness, Cd, Cf, Cm, curvature=curvature)
        next_free_state_index += dof(anext);
        push!(alast.children, anext); anext.parent = alast;
        alast = anext;
    end
    body = ArticulatedBody(a0, -1);

    dof(body, force=true);

    return body;
end




force_drag(fluid_density::Number, u::Vector{<:Number}) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> filament_drag(a.properties.ffi, i, fluid_density, u)
);
force_drag(fluid_density::Number, u::Function) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> filament_drag(a.properties.ffi, i, fluid_density, u(i.p, t))
);
function filament_drag(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number, u::Vector{<: Number})
    filament_normal_flow = antiproject(u .- i.V[4:6], i.X0[3,1:3]);
    
    return 0.5fluid_density * ffi.drag_coefficient * frontal_area(ffi, i, filament_normal_flow) * norm(filament_normal_flow) * filament_normal_flow;
end


force_skin_friction(fluid_density::Number, u::Vector{<:Number}) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> filament_skin_friction(a.properties.ffi, i, fluid_density, u)
);
force_skin_friction(fluid_density::Number, u::Function) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> filament_skin_friction(a.properties.ffi, i, fluid_density, u(i.p, t))
);
function filament_skin_friction(ffi::StandardFFIProperties, i::ArticulationHarness, fluid_density::Number, u::Vector{<: Number})
    filament_tangent_flow = project(u .- i.V[4:6], i.X0[3,1:3]);

    return 0.5fluid_density * ffi.friction_coefficient * skin_area(ffi, i, filament_tangent_flow) * norm(filament_tangent_flow) * filament_tangent_flow;
end


force_buoyancy(fluid_density::Number;g::Number=9.81, downwards::Vector{<:Number}=[0,0,-1]) = ExternalForce(
    (a::Articulation, i::ArticulationHarness, t::Real) -> -(fluid_density * a.properties.ffi.volume * g) .* downwards
);



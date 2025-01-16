
function inertia_matrix_about_origin(mass::Real, inertia_matrix_about_center::Matrix{<:Real}, OC::Vector{<:Real})
    return inertia_matrix_about_point(mass, inertia_matrix_about_center, OC, [0.0,0.0,0.0]);
end
function inertia_matrix_about_point(mass::Real, inertia_matrix_about_center::Matrix{<:Real}, OC::Vector{<:Real}, OP::Vector{<:Real})
    m = mass; PC = OC - OP;
    Ic = inertia_matrix_about_center;
    Itl = Ic + m*𝞦(PC)*𝞦(PC)';   
    Itr = m*𝞦(PC);
    Ibl = Itr';
    Ibr = diagm([m,m,m]);
    Io = [Itl Itr; Ibl Ibr];
    return Io;
end


CylinderArticulation(
    number::Integer, next_free_state_index::Integer, joint::JointType,
    parent::Articulation, length::Real, radius::Real, mass::Real, properties;
    curvature::Matrix{<:Real}=[1 0 0; 0 1 0; 0 0 1], pregeometry::Vector{<:Real}=[0,0,0]
) = Articulation(
    number, next_free_state_index, joint, parent, length,
    mass, cylinder_inertia(mass, length, radius), curvature, [0,0,length], pregeometry, [0,0,0.5length], properties
);

function cylinder_inertia(mass::Real, length::Real, radius::Real)
    inertia_xy = mass*(3*(radius^2) + (length^2))/12; inertia_z = mass * (radius^2) / 12;
    inertia = [inertia_xy 0 0; 0 inertia_xy 0; 0 0 inertia_z];
    return inertia;
end
function cylinder_mass_and_inertia(density::Real, length::Real, radius::Real)
    mass = density * length * pi * radius^2;
    return (mass, cylinder_inertia(mass, length, radius))
end
function rectangular_beam_inertia(mass::Real, length::Real, width::Real, thickness::Real)
    inertia_x = mass * ((length^2) + (width^2)) / 12;
    inertia_y = mass * ((length^2) + (thickness^2)) / 12;
    inertia_z = mass * ((thickness^2) + (width^2)) / 12;
    inertia = [inertia_x 0 0; 0 inertia_y 0; 0 0 inertia_z];
    return inertia;
end
function rectangular_beam_mass_and_inertia(density::Real, length::Real, width::Real, thickness::Real)
    mass = density * length * width * thickness;
    return (mass, rectangular_beam_inertia(mass, length, width, thickness))
end
function rectangular_beam_lumped_mass_and_inertia(density::Real, length::Real, width::Real, thickness::Real)
    mass = density * length * width * thickness;
    return (mass, lumped_mass_inertia(mass, thickness))
end
function lumped_mass_inertia(mass::Real, radius::Real)
    Ixyz = 0.4mass * (radius ^ 2);

    return [Ixyz 0 0; 0 Ixyz 0; 0 0 Ixyz];
end



function rotate_x(theta::Number)
    c = cos(theta); s = sin(theta);
    return [
        1   0   0;
        0  c  s;
        0 -s  c ]
end
function rotate_y(theta::Number)
    c = cos(theta); s = sin(theta);
    return [
        c   0 -s;
        0   1   0;
        s   0  c]
end
function rotate_z(theta::Number)
    c = cos(theta); s = sin(theta);
    return [
        c  s   0;
        -s  c   0;
        0   0   1]
end
function get_rotation(a::Vector{<: Number}, b::Vector{<: Number})
    # b = Ra
    # https://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d
    v = cross(a,b);
    c = dot(a,b);

    𝞦v = 𝞦(v);
    return I(3) + 𝞦v + 𝞦v*𝞦v/(1+c);
end


function get_rotation_b_is_ez(a::Vector{<: Number})
    v = [a[2], -a[1], 0];
    c = a[3];
    𝞦v = 𝞦(v);
    return I(3) + 𝞦v + 𝞦v*𝞦v/(1+c);
end
# Commonly used vector calculations

function project(unprojected::Vector{<:Number}, dir::Vector{<:Number})
    return (unprojected ⋅ dir) * dir;
end
function antiproject(unprojected::Vector{<:Number}, dir::Vector{<:Number})
    return unprojected - project(unprojected, dir);
end

# Use of ForwardDiff! flow_func must be auto-differentiable.
function flow_acceleration(flow::Function, x::Vector{<: Number}, t::Number)
    U_oft = (t_prime::Vector) -> flow(x, t_prime[1]);
    return vec(jacobian(U_oft, [t]))
end
function flow_jacobian(flow::Function, x::Vector{<: Number}, t::Number)
    U_ofx = (x_prime::Vector) -> flow(x_prime, t);
    return jacobian(U_ofx, x);
end
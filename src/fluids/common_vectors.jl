function project(unprojected::Vector{<:Number}, dir::Vector{<:Number})
    return (unprojected ⋅ dir) * dir;
end
function antiproject(unprojected::Vector{<:Number}, dir::Vector{<:Number})
    return unprojected - project(unprojected, dir);
end

function flow_acceleration(flow::Function, x::Vector{<: Number}, t::Number)
    U_oft = (t_prime::Vector) -> flow(x, t_prime[1]);
    return vec(jacobian(U_oft, [t]))
end
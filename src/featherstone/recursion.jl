
function forward_recurse(body::AbstractArticulatedBody, f::Function, initialstate, args...; retrace::Bool=false, kwargs...)
    return forward_recurse(body.articulation_zero.children[1], f, retrace, initialstate, args...; kwargs...);
end

function forward_recurse(a::Articulation, f::Function, retrace::Bool, state, args...; kwargs...)
    state = f(a, state, args...; kwargs...);
    for child in a.children
        state = forward_recurse(child, f, retrace, state, args...; kwargs...);
        if retrace
            state = f(a, state, args...; kwargs...);
        end
    end
    return state;
end

function forward_recurse!(body::AbstractArticulatedBody, f::Function, args...; retrace::Bool=false, kwargs...)
    forward_recurse!(body.articulation_zero.children[1], f, retrace, args...; kwargs...);
end

function forward_recurse!(a::Articulation, f::Function, retrace::Bool, args...; kwargs...)
    f(a, args...; kwargs...);
    for child in a.children
        forward_recurse!(child, f, retrace, args...; kwargs...)
        if retrace
            f(a, args...; kwargs...);
        end
    end
end


function forward_recurse_iterative!(body::AbstractArticulatedBody, f::Function, args...; kwargs...)

    to_visit::Array{Articulation} = [ body.articulation_zero.children[1] ];

    while !isempty(to_visit)
        next = pop!(to_visit);
        f(next, args...; kwargs...);
        append!(to_visit, next.children);
    end
end



function backward_recurse!(body::AbstractArticulatedBody, f::Function, args...; kwargs...)
    backward_recurse!(body.articulation_zero.children[1], f, args...; kwargs...);
end

function backward_recurse!(a::Articulation, f::Function, args...; kwargs...)
    for child in a.children
        backward_recurse!(child, f, args...; kwargs...)
    end
    f(a, args...; kwargs...);
end

function backward_recurse_iterative!(body::AbstractArticulatedBody, f::Function, args...; kwargs...)
    to_visit::Array{Articulation} = [ body.articulation_zero.children[1] ];
    to_apply::Array{Articulation} = [];

    while !isempty(to_visit)
        next = pop!(to_visit);
        push!(to_apply, next);
        append!(to_visit, next.children);
    end

    while !isempty(to_apply)
        next = pop!(to_apply);
        f(next, args...; kwargs...);
    end
end

function forward_recurse(body::AbstractArticulatedBody, f::Function, initialstate, args...; retrace::Bool=false, kwargs...)
    return forward_recurse(body.articulation_zero.children[1], f, retrace, initialstate, args...; kwargs...);
end

function forward_recurse(a::Articulation, f::Function, retrace::Bool, state, args...; kwargs...)
    state = f(a, state, args...; kwargs...);
    for child in a.children
        state = forward_recurse(child, f, retrace, state, args...; kwargs...);
    end
    if retrace
        state = f(a, state, args...; kwargs...);
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
    end
    if retrace
        f(a, args...; kwargs...);
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

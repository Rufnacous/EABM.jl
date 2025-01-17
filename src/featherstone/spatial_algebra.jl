# See Chapter 2 of Rigid Body Dynamics Algorithms by Roy Featherstone (2008).

function 𝞦(x::Vector{<: Real})
    if length(x) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; x is length %d.", length(x))));
    end
    return [ 0    -x[3]  x[2] ;
             x[3]  0    -x[1] ;
            -x[2]  x[1]  0    ];
end

# Spatial cross product
function ⨱(a::Vector{<: Real}, b::Vector{<: Real})
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    return [ 𝞦(a[1:3]) zeros(3,3) ; 𝞦(a[4:6]) 𝞦(a[1:3]) ] * b;
end
# Dual of spatial cross product
function ⨳(a::Vector{<: Real}, b::Vector{<: Real})
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    return [ 𝞦(a[1:3]) 𝞦(a[4:6]) ; zeros(3,3) 𝞦(a[1:3]) ] * b;
end
function xlt(r::Vector{<:Real})
    if length(r) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; r is length %d.", length(r))));
    end
    XT = Matrix{typeof(r).parameters[1]}(diagm(ones(6)))
    XT[4:6,1:3] = -𝞦(r);
    return XT;
end

function diag_spatial(m::Matrix{<:Real})
    if ndims(m) != 2
        throw(DimensionMismatch(@sprintf("diag_spatial is only defined for matrices of size 3x3 and m has %d dims.", ndims(m))));
    end
    if size(m) != (3,3)
        throw(DimensionMismatch(@sprintf("diag_spatial is only defined for matrices of size 3x3 and m is of size %dx%d.",size(m)...)));
    end
    return [m zeros(3,3); zeros(3,3) m];
end
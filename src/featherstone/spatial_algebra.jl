# See Chapter 2 of Rigid Body Dynamics Algorithms by Roy Featherstone (2008).

function 𝞦(x::Vector{<: Real})
    if length(x) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; x is length %d.", length(x))));
    end
    return [ 0    -x[3]  x[2] ;
             x[3]  0    -x[1] ;
            -x[2]  x[1]  0    ];
end
function 𝞦!(𝞦p,x::AbstractVector{<: Real})
    if length(x) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; x is length %d.", length(x))));
    end
    𝞦p[1,1] = 0; 𝞦p[2,2] = 0; 𝞦p[3,3] = 0;
    𝞦p[1,2] = -x[3]; 𝞦p[1,3] = x[2];
    𝞦p[2,1] = x[3]; 𝞦p[2,3] = -x[1];
    𝞦p[3,1] = -x[2]; 𝞦p[3,2] = x[1];
end
function 𝞦_minus!(𝞦p,x::AbstractVector{<: Real})
    if length(x) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; x is length %d.", length(x))));
    end
    𝞦p[1,1] = 0; 𝞦p[2,2] = 0; 𝞦p[3,3] = 0;
    𝞦p[1,2] = x[3]; 𝞦p[1,3] = -x[2];
    𝞦p[2,1] = -x[3]; 𝞦p[2,3] = x[1];
    𝞦p[3,1] = x[2]; 𝞦p[3,2] = -x[1];
end

function cross!(c,a,b)
    c[1] = a[2]*b[3] - a[3]*b[2];
    c[2] = a[3]*b[1] - a[1]*b[3];
    c[3] = a[1]*b[2] - a[2]*b[1];
end
# Spatial cross product
function spatial_cross(a,b)
    return a ⨱ b;
end
function ⨱(a::Vector{<: Real}, b::Vector{<: Real})
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    return [ 𝞦(a[1:3]) zeros(3,3) ; 𝞦(a[4:6]) 𝞦(a[1:3]) ] * b;
end
function spatial_cross!(c, a, b)
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end

    a_top = @view a[1:3];
    a_btm = @view a[4:6];
    b_top = @view b[1:3];
    b_btm = @view b[4:6];
    c_top = @view c[1:3];
    c_btm = @view c[4:6];

    cross!(c_btm, a_btm, b_top);
    cross!(c_top, a_top, b_btm);
    c_btm .+= c_top;

    cross!(c_top, a_top, b_top);
end

# Dual of spatial cross product
function dual_spatial_cross(a,b)
    return a ⨳ b;
end
function ⨳(a::Vector{<: Real}, b::Vector{<: Real})
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    return [ 𝞦(a[1:3]) 𝞦(a[4:6]) ; zeros(3,3) 𝞦(a[1:3]) ] * b;
end
function dual_spatial_cross!(c, a, b)
    if length(a) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end
    if length(b) != 6
        throw(DimensionMismatch(@sprintf("⨱ is only defined for vectors of length 6; a and b are lengths %d and %d.",length(a), length(b))));
    end

    a_top = @view a[1:3];
    a_btm = @view a[4:6];
    b_top = @view b[1:3];
    b_btm = @view b[4:6];
    c_top = @view c[1:3];
    c_btm = @view c[4:6];

    cross!(c_top, a_top, b_top);
    cross!(c_btm, a_btm, b_btm);
    c_top .+= c_btm;

    cross!(c_btm, a_top, b_btm);
end
function xlt(r::Vector{<:Real})
    if length(r) != 3
        throw(DimensionMismatch(@sprintf("𝞦 is only defined for vectors of length 3; r is length %d.", length(r))));
    end
    XT = Matrix{typeof(r).parameters[1]}(diagm(ones(6)))
    XT[4:6,1:3] = -𝞦(r);
    return XT;
end
function xlt!(A::AbstractMatrix{<:Real}, r::Vector{<:Real})
    A .= 0;
    A[1] = 1;
    A[2] = 1;
    A[3] = 1;
    A[4] = 1;
    A[5] = 1;
    A[6] = 1;
    Abl = @view A[4:6,1:3];
    𝞦_minus!(Abl, r);
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

function xltT_mul_btm!(pA, com, X0, fx_o, 𝞦com)

    pAtop = @view pA[1:3];
    pAbtm = @view pA[4:6];

    X0tr = @view X0[1:3,4:6];
    𝞦!(𝞦com, com);
    mul!(pAtop, (X0tr + 𝞦com), fx_o);

    X0br = @view X0[4:6,4:6];
    mul!(pAbtm, X0br,fx_o);

end
"""
Rotation matrix for a 2D truss element

Cx, Cy are the directional cosines of the element in global coordinates

Cx, Cy = (endPosition - startPosition) / length

IE they are the components of the normalized vector that represents the element.
"""
function R2DTruss(Cx::Real, Cy::Real)

    @assert norm([Cx, Cy]) ≈ 1

    return [Cx Cy 0 0 ; 0 0 Cx Cy]
end

"""
Rotation matrix for a 2D truss using start and end points

positionStart = [xstart, ystart]: position of starting node
positionEnd = [xend, yend]: position of ending node
L: length of element
"""
function R2DTruss(positionStart::Vector{Real}, positionEnd::Vector{Real})

    Cx, Cy = normalize(positionEnd .- positionStart)

    return [Cx Cy 0 0 ; 0 0 Cx Cy]
end

"""
Rotation matrix for a 3D truss element
Cx, Cy, Cz are the directional cosines of the element in global coordinates
Cx, Cy, Cz = (endPosition - startPosition) / length
"""
function R3DTruss(Cx::Real, Cy::Real, Cz::Real)

    @assert norm([Cx, Cy, Cz]) ≈ 1

    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]

end

"""
Rotation matrix for a 3D truss using start and end points

positionStart = [xstart, ystart, zstart]: position of starting node
positionEnd = [xend, yend, zend]: position of ending node
L: length of element
"""
function R3DTruss(positionStart::Vector{Real}, positionEnd::Vector{Real})

    Cx, Cy, Cz = normalize(positionEnd .- positionStart)

    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

"""
Rotation matrix for a 2D frame element

Cx, Cy are the directional cosines of the element in global coordinates

Cx, Cy = (endPosition - startPosition) / length

IE they are the components of the normalized vector that represents the element.
"""
function R2DFrame(Cx::Real, Cy::Real)
    
    @assert norm([Cx, Cy]) ≈ 1

    R = [Cx Cy 0 0 0 0;
        -Cy Cx 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 Cx Cy 0;
        0 0 0 -Cy Cx 0;
        0 0 0 0 0 1]
    
    return R
end

"""
Rotation matrix for a 2D frame element using start and end points

positionStart = [xstart, ystart]: position of starting node
positionEnd = [xend, yend]: position of ending node
L: length of element
"""
function R2DFrame(positionStart::Vector{Real}, positionEnd::Vector{Real})
    
    Cx, Cy = normalize(positionEnd .- positionStart)

    R = [Cx Cy 0 0 0 0;
        -Cy Cx 0 0 0 0;
        0 0 1 0 0 0;
        0 0 0 Cx Cy 0;
        0 0 0 -Cy Cx 0;
        0 0 0 0 0 1]
    
    return R
end

"""
Rotation matrix for a 3D frame element
Cx, Cy, Cz are the directional cosines of the element in global coordinates
Cx, Cy, Cz = (endPosition - startPosition) / length

Ψ is the roll angle of the element with respect to local x axis
By default, Ψ = π/2, which ensures that the local Z axis is parallel to the XY ground plane, IE the strong axis is fully utilized under global gravity loading
"""
function R3DFrame(Cx::Real, Cy::Real, Cz::Real; Ψ = pi/2)

    xLocal = [Cx, Cy, Cz]

    @assert norm(xLocal) ≈ 1

    if norm(cross(xLocal, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. Cy 0.;
            -Cy*cos(Ψ) 0 sin(Ψ);
            Cy*sin(Ψ) 0 cos(Ψ)]
    else # all other
        b1 = (-Cx * Cy * cos(Ψ) - Cz * sin(Ψ)) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cos(Ψ)
        b3 = (-Cy * Cz * cos(Ψ) + Cx * sin(Ψ)) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sin(Ψ) - Cz * cos(Ψ)) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sin(Ψ)
        c3 = (Cy * Cz * sin(Ψ) + Cx * cos(Ψ)) / sqrt(Cx^2 + Cz^2)

        Λ = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end

"""
Rotation matrix for a 3D frame element
positionStart = [xstart, ystart, zstart]: position of starting node
positionEnd = [xend, yend, zend]: position of ending node
L: length of element

Ψ is the roll angle of the element with respect to local x axis
By default, Ψ = π/2, which ensures that the local Z axis is parallel to the XY ground plane, IE the strong axis is fully utilized under global gravity loading
"""
function R3DFrame(positionStart::Vector{Real}, positionEnd::Vector{Real}; Ψ = pi/2)
    
    xLocal = normalize(positionEnd .- positionStart)
    Cx, Cy, Cz = xLocal

    if norm(cross(xLocal, globalY)) < tol #special case for horizontal members aligned with global Y
        Λ = [0. Cy 0.;
            -Cy*cos(Ψ) 0 sin(Ψ);
            Cy*sin(Ψ) 0 cos(Ψ)]
    else # all other
        b1 = (-Cx * Cy * cos(Ψ) - Cz * sin(Ψ)) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cos(Ψ)
        b3 = (-Cy * Cz * cos(Ψ) + Cx * sin(Ψ)) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sin(Ψ) - Cz * cos(Ψ)) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sin(Ψ)
        c3 = (Cy * Cz * sin(Ψ) + Cx * cos(Ψ)) / sqrt(Cx^2 + Cz^2)

        Λ = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [Λ zeros(3,9); zeros(3,3) Λ zeros(3,6); zeros(3,6) Λ zeros(3,3); zeros(3,9) Λ]

    return R
end
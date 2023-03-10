"""
Rotation matrix for a 2D truss element

Cx, Cy are the directional cosines of the element in global coordinates

Cx, Cy = (endPosition - startPosition) / length

IE they are the components of the normalized vector that represents the element.
"""
function R2DTruss(Cx::Real, Cy::Real)
    return [Cx Cy 0 0 ; 0 0 Cx Cy]
end

"""
Rotation matrix for a 2D truss using start and end points

positionStart = [xstart, ystart]: position of starting node
positionEnd = [xend, yend]: position of ending node
L: length of element
"""
function R2DTruss(positionStart::Vector{Real}, positionEnd::Vector{Real}, L::Real)

    Cx, Cy = (positionEnd .- positionStart) ./ L

    return [Cx Cy 0 0 ; 0 0 Cx Cy]
end

"""
Rotation matrix for a 3D truss element
Cx, Cy, Cz are the directional cosines of the element in global coordinates
Cx, Cy, Cz = (endPosition - startPosition) / length
"""
function R3DTruss(Cx::Real, Cy::Real, Cz::Real)

    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]

end

"""
Rotation matrix for a 3D truss using start and end points

positionStart = [xstart, ystart, zstart]: position of starting node
positionEnd = [xend, yend, zend]: position of ending node
L: length of element
"""
function R2DTruss(positionStart::Vector{Real}, positionEnd::Vector{Real}, L::Real)

    Cx, Cy, Cz = (positionEnd .- positionStart) ./ L

    return [Cx Cy Cz 0 0 0; 0 0 0 Cx Cy Cz]
end

"""
Rotation matrix for a 2D frame element

Cx, Cy are the directional cosines of the element in global coordinates

Cx, Cy = (endPosition - startPosition) / length

IE they are the components of the normalized vector that represents the element.
"""
function R2DFrame(Cx::Real, Cy::Real)
    
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
function R2DFrame(positionStart::Vector{Real}, positionEnd::Vector{Real}, L::Real)
    
    Cx, Cy = (positionEnd .- positionStart) ./ L

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

?? is the roll angle of the element with respect to local x axis
By default, ?? = ??/2, which ensures that the local Z axis is parallel to the XY ground plane, IE the strong axis is fully utilized under global gravity loading
"""
function R3DFrame(Cx::Real, Cy::Real, Cz::Real; ?? = pi/2)

    xLocal = [Cx, Cy, Cz]

    if norm(cross(xLocal, globalY)) < tol #special case for horizontal members aligned with global Y
        ?? = [0. Cy 0.;
            -Cy*cos(??) 0 sin(??);
            Cy*sin(??) 0 cos(??)]
    else # all other
        b1 = (-Cx * Cy * cos(??) - Cz * sin(??)) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cos(??)
        b3 = (-Cy * Cz * cos(??) + Cx * sin(??)) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sin(??) - Cz * cos(??)) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sin(??)
        c3 = (Cy * Cz * sin(??) + Cx * cos(??)) / sqrt(Cx^2 + Cz^2)

        ?? = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [?? zeros(3,9); zeros(3,3) ?? zeros(3,6); zeros(3,6) ?? zeros(3,3); zeros(3,9) ??]

    return R
end

"""
Rotation matrix for a 3D frame element
positionStart = [xstart, ystart, zstart]: position of starting node
positionEnd = [xend, yend, zend]: position of ending node
L: length of element

?? is the roll angle of the element with respect to local x axis
By default, ?? = ??/2, which ensures that the local Z axis is parallel to the XY ground plane, IE the strong axis is fully utilized under global gravity loading
"""
function R3DFrame(positionStart::Vector{Real}, positionEnd::Vector{Real}, L::Real; ?? = pi/2)
    
    xLocal = normalize(positionEnd .- positionStart)
    Cx, Cy, Cz = xLocal

    if norm(cross(xLocal, globalY)) < tol #special case for horizontal members aligned with global Y
        ?? = [0. Cy 0.;
            -Cy*cos(??) 0 sin(??);
            Cy*sin(??) 0 cos(??)]
    else # all other
        b1 = (-Cx * Cy * cos(??) - Cz * sin(??)) / sqrt(Cx^2 + Cz^2)
        b2 = sqrt(Cx^2 + Cz^2) * cos(??)
        b3 = (-Cy * Cz * cos(??) + Cx * sin(??)) / sqrt(Cx^2 + Cz^2)

        c1 = (Cx * Cy * sin(??) - Cz * cos(??)) / sqrt(Cx^2 + Cz^2)
        c2 = -sqrt(Cx^2 + Cz^2) * sin(??)
        c3 = (Cy * Cz * sin(??) + Cx * cos(??)) / sqrt(Cx^2 + Cz^2)

        ?? = [Cx Cy Cz; 
            b1 b2 b3; 
            c1 c2 c3]
    end
    
    R = [?? zeros(3,9); zeros(3,3) ?? zeros(3,6); zeros(3,6) ?? zeros(3,3); zeros(3,9) ??]

    return R
end
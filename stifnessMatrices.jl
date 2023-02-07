
"""
Truss element

E: Modulus of Elasticity
A: Cross section area
L: Element length
"""
function KTruss(E::Real, A::Real, L::Real)
    return E * A / L .* [1 -1; -1 1]
end

"""
3D beam element

E: Modulus of Elasticity
A: Cross section area
L: Element length
G: Shear Modulus
Izz: Moment of Inertia in local Z (strong)
Iyy: Moment of Inertia in local Y (weak)
J: Torsional Constant
"""
function K3DFrame(E::Real, A::Real, L::Real, G::Real, Izz::Real, Iyy::Real, J::Real)
    k = E / L^3 * [
        A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 12Izz 0 0 0 6L*Izz 0 -12Izz 0 0 0 6L*Izz;
        0 0 12Iyy 0 -6L*Iyy 0 0 0 -12Iyy 0 -6L*Iyy 0;
        0 0 0 G*J*L^2/E 0 0 0 0 0 -G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 4L^2*Iyy 0 0 0 6L*Iyy 0 2L^2*Iyy 0;
        0 6L*Izz 0 0 0 4L^2*Izz 0 -6L*Izz 0 0 0 2L^2*Izz;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -12Izz 0 0 0 -6L*Izz 0 12Izz 0 0 0 -6L*Izz;
        0 0 -12Iyy 0 6L*Iyy 0 0 0 12Iyy 0 6L*Iyy 0;
        0 0 0 -G*J*L^2/E 0 0 0 0 0 G*J*L^2/E 0 0;
        0 0 -6L*Iyy 0 2L^2*Iyy 0 0 0 6L*Iyy 0 4L^2*Iyy 0;
        0 6L*Izz 0 0 0 2L^2*Izz 0 -6L*Izz 0 0 0 4L^2*Izz
        ]

    return k
end

"""
3D beam element with a hinge at the starting node

E: Modulus of Elasticity
A: Cross section area
L: Element length
Izz: Moment of Inertia in local Z (strong)
Iyy: Moment of Inertia in local Y (weak)
"""
function K3DFrame_FreeFixed(E::Real, A::Real, L::Real, Izz::Real, Iyy::Real)

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Izz 0 0 0 0 0 -3Izz 0 0 0 3L*Izz;
        0 0 3Iyy 0 0 0 0 0 -3Iyy 0 -3L*Iyy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Izz 0 0 0 0 0 3Izz 0 0 0 -3L*Izz;
        0 0 -3Iyy 0 0 0 0 0 3Iyy 0 3L*Iyy 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iyy 0 0 0 0 0 3L*Iyy 0 3L^2*Iyy 0;
        0 3L*Izz 0 0 0 0 0 -3L*Izz 0 0 0 3L^2*Izz    
    ]

    return k
end

"""
3D beam element with a hinge at the ending node

E: Modulus of Elasticity
A: Cross section area
L: Element length
Izz: Moment of Inertia in local Z (strong)
Iyy: Moment of Inertia in local Y (weak)
"""
function K3DFrame_FixedFree(E::Real, A::Real, L::Real, Izz::Real, Iyy::Real)

    k = E / L^3 .* [A*L^2 0 0 0 0 0 -A*L^2 0 0 0 0 0;
        0 3Izz 0 0 0 3L*Izz 0 -3Izz 0 0 0 0;
        0 0 3Iyy 0 -3L*Iyy 0 0 0 -3Iyy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 -3L*Iyy 0 3L^2*Iyy 0 0 0 3L*Iyy 0 0 0;
        0 3L*Izz 0 0 0 3L^2*Izz 0 -3L*Izz 0 0 0 0 ;
        -A*L^2 0 0 0 0 0 A*L^2 0 0 0 0 0;
        0 -3Izz 0 0 0 -3L*Izz 0 3Izz 0 0 0 0;
        0 0 -3Iyy 0 3L*Iyy 0 0 0 3Iyy 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 0 0    
    ]

    return k
end

"""
3D beam element with a hinges at both ends (truss element) -- use sparingly

E: Modulus of Elasticity
A: Cross section area
L: Element length
"""
function K3DFrame_FreeFree(E::Real, A::Real, L::Real)

    k = zeros(12, 12)
    k[1, 1] = k[7, 7] = 1
    k[1, 7] = k[7, 1] = -1

    return E * A / L .* k
end

"""
2D beam element

E: Modulus of Elasticity
A: Cross section area
L: Element length
I: Moment of Inertia in out-of-plane axis
"""
function K2DFrame(E::Real, A::Real, L::Real, I::Real)

    k = E * I / L^3 .* [A*L^2/I 0 0 -A*L^2/I 0 0;
        0 12 6L 0 -12 6L;
        0 6L 4L^2 0 -6L 2L^2;
        -A*L^2/I 0 0 A*L^2/I 0 0;
        0 -12 -6L 0 12 -6L;
        0 6L 2L^2 0 -6L 4L^2]

    return k
end

"""
2D beam element with hinge at starting node

E: Modulus of Elasticity
A: Cross section area
L: Element length
I: Moment of Inertia in out-of-plane axis
"""
function K2DFrame_FreeFixed(E::Real, A::Real, L::Real, I::Real)

    k = E * I / L^3 .* [A*L^2/I 0 0 -A*L^2/I 0 0;
        0 3 0 0 -3 3L;
        0 0 0 0 0 0;
        -A*L^2/I 0 0 A*L^2/I 0 0;
        0 -3 0 0 3 -3L;
        0 3L 0 0 -3L 3L^2]

    return k
end

"""
2D beam element with hinge at ending node

E: Modulus of Elasticity
A: Cross section area
L: Element length
I: Moment of Inertia in out-of-plane axis
"""
function K2DFrame_FixedFree(E::Real, A::Real, L::Real, I::Real)

    k = E * I / L^3 .* [A*L^2/I 0 0 -A*L^2/I 0 0;
        0 3 3L 0 -3 0;
        0 3L 3L^2 0 -3L 0;
        -A*L^2/I 0 0 A*L^2/I 0 0;
        0 -3 -3L 0 3 0;
        0 0 0 0 0 0]

    return k
end

"""
2D beam element with hinges at both ends (ie a truss)

E: Modulus of Elasticity
A: Cross section area
L: Element length
I: Moment of Inertia in out-of-plane axis
"""
function K2DFrame_FixedFree(E::Real, A::Real, L::Real)

    k = zeros(6, 6)
    k[1, 1] = k[4, 4] = 1
    k[1, 4] = k[4, 1] = -1

    return E * A / L .* k
end
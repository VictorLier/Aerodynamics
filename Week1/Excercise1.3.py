import numpy as np
import sympy as sp

Tau, a = sp.symbols('Tau, a')


def VelocityInduction(R1,R2):
    R0 = R1 - R2
    r1 = sp.sqrt(R1[0]**2 + R1[1]**2)
    r2 = sp.sqrt(R2[0]**2 + R2[1]**2)
    R0ting = R0.T*(R1/r1-R2/r2)
    Q = Tau / (4*sp.pi) * R1.cross(R2) / ((R1.cross(R2)).norm()**2) *R0ting
    return Q

V1 = sp.Matrix([-a/2,-a/2,0])
V2 = sp.Matrix([-a/2,a/2,0])
V3 = sp.Matrix([a/2,a/2,0])
V4 = sp.Matrix([a/2,-a/2,0])

Q12 = VelocityInduction(V1,V2)

print(Q12)

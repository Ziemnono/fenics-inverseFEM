#!/usr/bin/env python
#description     :Perform an incompressible hyper-elastic simulation, in this
#				  programm, we reproduce the generalized shear experiment
#authors         :Arnaud Mazier & Jack Hale
#contact         :mazier.arnaud@gmail.com
#date            :10/2020
#python_version  :3.8.2
#paper			 :
#==============================================================================

import numpy as np
from dolfin import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"

mesh_size = 5
mesh = UnitCubeMesh(mesh_size, mesh_size, mesh_size)
with XDMFFile('results/initial_configuration_generalized_shear.xdmf') as xdmf:
    xdmf.write(mesh)

V = VectorFunctionSpace(mesh, "CG", 1)

C_01 = 10
C_10 = 1
D = 1000000000.0
k = 1.0

u_ = Function(V)
u = TrialFunction(V)
v = TestFunction(V)

# kinematics with deformation gradient, right Green-Cauchy and Lagrangian strain tensor
I = Identity(3)
F = variable(grad(u_) + I)
C = variable(F.T * F)
E = variable(0.5 * (C - I))

# invariants and Jacobian
I_C = tr(C)
II_C = (1.0 / 2.0) * (tr(C) ** 2 - tr(C * C))
III_C = det(C)
J = III_C ** (1.0 / 2.0)

# stored strain energy density (incompressible Mooney-Rivlin model)
psi = Constant(C_01) * (J ** (-2.0 / 3.0) * I_C - 3) + Constant(C_10) * (J ** (-4.0 / 3.0) * II_C - 3) + Constant(D) * (
        J - 1) ** 2

# second Piola-Kirchoff conjugate with incremental Green-Lagrange
S = 2.0 * diff(psi, C)
D_E = derivative(E, u_, v)
G = inner(S, D_E) * dx
J = derivative(G, u_, u)


def top(x, on_boundary):
    return near(x[1], 1.0) and on_boundary


def bottom(x, on_boundary):
    return near(x[1], 0.0) and on_boundary


def right(x, on_boundary):
    return near(x[0], 1.0) and on_boundary


def left(x, on_boundary):
    return near(x[0], 0.0) and on_boundary


r = Expression(("k*x[1]*x[1]", "0.0", "0.0"), k=k, degree=1)
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), bottom)] + \
      [DirichletBC(V, r, top)] + \
      [DirichletBC(V, r, right)] + \
      [DirichletBC(V, r, left)]

bc_index = []
for vertex in vertices(mesh):
    if vertex.point().x() == 0.0 or vertex.point().x() == 1.0:
        bc_index.append(vertex.index())
np.savetxt("results/generalized_shear_bc_index.txt", bc_index, delimiter=",")

solver_parameters = {"newton_solver": {"linear_solver": "mumps"}}

solve(G == 0, u_, J=J, bcs=bcs, solver_parameters=solver_parameters)
ALE.move(mesh, u_)
with XDMFFile("results/forward_configuration_generalized_shear.xdmf") as xdmf:
    xdmf.write(mesh)

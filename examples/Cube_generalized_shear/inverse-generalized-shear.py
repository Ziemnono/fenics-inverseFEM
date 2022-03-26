#!/usr/bin/env python
# description     :Perform an inverse incompressible hyper-elastic simulation,
#				  in this programm, we apply an inverse displacement on a 3D
#				  cube previously deformed to recover its initial shape
# authors         :Arnaud Mazier & Jack Hale
# contact         :mazier.arnaud@gmail.com
# date            :10/2020
# python_version  :3.8.2
# paper			 :
# ==============================================================================

import numpy as np
from dolfin import *

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"

mesh = Mesh()
with XDMFFile("results/forward_configuration_generalized_shear.xdmf") as infile:
    infile.read(mesh)

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
F = variable(inv(grad(u_) + I))
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


r = Expression(("-k*x[1]*x[1]", "0.0", "0.0"), k=k, degree=1)
bcs = [DirichletBC(V, Constant((0.0, 0.0, 0.0)), bottom)] + \
      [DirichletBC(V, r, top)]


def dirichlet_bc_on_vertices(mesh, space, value, index_number):
    string_list = ""
    for vertex in vertices(mesh):
        if vertex.index() == index_number:
            string_list = string_list + "near(x[0]," + \
                          str(vertex.point().x()) + \
                          ") && near(x[1]," \
                          + str(vertex.point().y()) + \
                          ") && near(x[2]," + \
                          str(vertex.point().z()) + ") || "
    return DirichletBC(space, value, string_list[0:-3], method="pointwise")


indices = np.loadtxt("results/generalized_shear_bc_index.txt").astype(int).tolist()
for i in indices:
    bcs.append(dirichlet_bc_on_vertices(mesh, V, r, i))

solver_parameters = {"newton_solver": {"linear_solver": "mumps"}}

solve(G == 0, u_, J=J, bcs=bcs, solver_parameters=solver_parameters)
ALE.move(mesh, u_)
with XDMFFile("results/inverse_configuration_generalized_shear.xdmf") as xdmf:
    xdmf.write(mesh)

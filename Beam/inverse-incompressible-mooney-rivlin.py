#!/usr/bin/env python
# description     :Perform an inverse incompressible hyper-elastic simulation,
#				  in this programm, we apply the inverse gravity on a soft beam
#				  previously deformed with gravity to recover its initial shape
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
parameters["form_compiler"]["quadrature_degree"] = 2

# loading the defomed geometry
mesh_file = "Beam2k"
mesh = Mesh()
with XDMFFile("results/forward_deformed_configuration" + mesh_file + ".xdmf") as infile:
    infile.read(mesh)

# mixed function space, 2nd order for displacement and 1st order for pressure
V = VectorElement("CG", mesh.ufl_cell(), 2)
P = FiniteElement("CG", mesh.ufl_cell(), 1)
U = FunctionSpace(mesh, MixedElement([V, P]))

# material properties
C_01 = 101709.668
C_10 = 151065.460
D = 7.965272689E08
rho = 963.0
rho_e = Expression("rho", rho=rho, degree=0)

# loading, gravity in this case
g = 9.81

u_p_ = Function(U)
u = TrialFunction(U)
v = TestFunction(U)

u_t, p_t = split(v)
u_, p_ = split(u_p_)

# kinematics with inversed deformation gradient, right Green-Cauchy and 
# Lagrangian strain tensor 
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
psi = Constant(C_01) * (J ** (-2.0 / 3.0) * I_C - 3) + Constant(C_10) * (J ** (-4.0 / 3.0) * II_C - 3) + p_ * (
            J - 1) - (1.0 / (4.0 * Constant(D))) * p_ ** 2

# second Piola-Kirchoff conjugate with incremental Green-Lagrange
S = 2.0 * diff(psi, C)
sigma = J ** (-1.0) * F * S * F.T
D_E = derivative(E, u_p_, v)
G = inner(sigma, sym(grad(u_t))) * dx + inner(((J - 1) - p_ / Constant(2.0 * D)), p_t) * dx - inner(Constant(g) * rho_e,
                                                                                                    -u_t[1]) * dx
J = derivative(G, u_p_, u)


# define the ROI (Region Of Interest) corresponding to the left side of the beam
def left(x, on_boundary):
    return near(x[0], 0.0) and on_boundary


# apply null Dirichlet boundary conditions on the ROI
bcs = [DirichletBC(U.sub(0), Constant((0.0, 0.0, 0.0)), left)]

solver_parameters = {"newton_solver": {"linear_solver": "mumps"}}

# for incremental load (slower than direct but improve convergence chances)
import numpy as np

steps = 10
rhos = np.linspace(0.0, rho, num=steps)
for rho in rhos:
    rho_e.rho = rho
    solve(G == 0, u_p_, J=J, bcs=bcs, solver_parameters=solver_parameters)

## direct load (quicker but has more chances to diverge)
# solve(G == 0, u_p_, J=J, bcs=bcs,solver_parameters = solver_parameters)

# export the deformed mesh in the results directory
ALE.move(mesh, u_p_.sub(0))
with XDMFFile("results/inverse_configuration" + mesh_file + ".xdmf") as xdmf:
    xdmf.write(mesh)

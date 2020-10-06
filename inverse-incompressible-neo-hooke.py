#!/usr/bin/env python
#description     :Perform an inverse incompressible hyper-elastic simulation,
#				  in this programm, we apply an inverse displacement on a 2D beam
#				  square previously deformed to recover its initial shape
#authors         :Arnaud Mazier & Jack Hale
#contact         :mazier.arnaud@gmail.com
#date            :10/2020
#python_version  :3.8.2
#paper			 : 
#==============================================================================

from dolfin import *
import numpy as np

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"
parameters["form_compiler"]["quadrature_degree"] = 2


# loading the deformed geometry
mesh = Mesh()
with XDMFFile("results/forward_deformed_configuration_cube.xdmf") as infile:
    infile.read(mesh)

# mixed function space, 2nd order for displacement and 1st order for pressure
V = VectorElement("CG", mesh.ufl_cell(), 2) 
P = FiniteElement("CG", mesh.ufl_cell(), 1)
U = FunctionSpace(mesh, MixedElement([V,P])) 

# material properties
mu = 0.6
lmbda = 10000.0 

# loading, gravity in this case
g = 9.81

u_p_ = Function(U)
u = TrialFunction(U)
v = TestFunction(U)

u_t, p_t = split(v)
u_, p_ = split(u_p_) 

# kinematics with deformation gradient, right Green-Cauchy and Lagrangian strain tensor 
I = Identity(2)
F = variable(inv(grad(u_) + I))
C = variable(F.T*F)
E = variable(0.5*(C - I))

# invariariants and Jacobian
I_C = tr(C)
II_C = (1.0/2.0)*(tr(C)**2 - tr(C*C))
III_C = det(C)
J = III_C**(1.0/2.0)

# stored strain energy density (incompressible Neo-Hooke model)
psi = (Constant(mu)/2.0)*(I_C - 2) - mu*ln(J) + p_*ln(J) - (1.0/(2.0*Constant(lmbda)))*p_**2

# second Piola-Kirchoff conjugate with incremental Green-Lagrange
S = 2.0*diff(psi, C)
D_E = derivative(E, u_p_, v)
G = inner(S, D_E)*dx + inner((ln(J) - p_/Constant(lmbda)), p_t)*dx
J = derivative(G, u_p_, u)

# define the ROI (Region Of Interest) corresponding to the top and bottom of the square 
def top(x, on_boundary):
     return near(x[1], 0.8) and on_boundary

def bottom(x, on_boundary):
     return near(x[1], -1.0) and on_boundary

bcs = [DirichletBC(U.sub(0), Constant((0.0, 0.0)), bottom)] + \
      [DirichletBC(U.sub(0), Constant((0.0, 0.2)), top)]


# direct solver using snes solver
solver_parameters = {"nonlinear_solver": "snes",
                     "snes_solver":
                     {"linear_solver": "lu",
                     "line_search": "bt",
                     "maximum_iterations": 1000}}
solve(G == 0, u_p_, J=J, bcs=bcs,solver_parameters = solver_parameters)

# export the deformed mesh in the results directory
ALE.move(mesh, u_p_.sub(0))
with XDMFFile("results/inverse_initial_configuration_cube.xdmf") as xdmf:
    xdmf.write(mesh)

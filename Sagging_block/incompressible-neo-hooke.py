#!/usr/bin/env python
# description     :Perform Monte Carlo simulation for a 3D sagging cube, in this
#                 programm, we run multiple simulations with random parameters 
#                 for a forward and inverse analysis to recover its initial shape. 
# authors         :Arnaud Mazier & Jack Hale
# contact         :mazier.arnaud@gmail.com
# date            :10/2021
# python_version  :3.8.10
# paper          :
# ==============================================================================

import numpy as np
from dolfin import *
import sys

parameter1 =  np.random.uniform(2.0e5, 2.0e3, 100)
parameter2 =  np.random.uniform(8.0e6, 8.0e4, 100)

def main():

    param_1 = []
    param_2 = []
    total_mse = []
    for i in parameter1 :
        for j in parameter2:
            print("=========", i, "    " , j, "=======")

            ### FORWARD  ###

            parameters["form_compiler"]["cpp_optimize"] = True
            parameters["form_compiler"]["representation"] = "uflacs"
            parameters["form_compiler"]["cpp_optimize_flags"] = "-O3 -ffast-math -march=native"
            parameters["form_compiler"]["quadrature_degree"] = 2

            # creating a 3D cube geometry
            mesh = UnitCubeMesh(5, 5, 5)
            mesh.coordinates()[:, :] *= 2.0
            mesh.coordinates()[:, 0] -= 1.0
            mesh.coordinates()[:, 1] -= 1.0
            mesh.coordinates()[:, 2] -= 1.0

            # save the geometry in the result directory
            with XDMFFile("results/initial_configuration_cube.xdmf") as xdmf:
                xdmf.write(mesh)

            # mixed function space, 2nd order for displacement and 1st order for pressure
            V = VectorElement("CG", mesh.ufl_cell(), 2)
            P = FiniteElement("CG", mesh.ufl_cell(), 1)
            U = FunctionSpace(mesh, MixedElement([V, P]))

            # material properties
            mu = Constant(i)
            lmbda = Constant(j)
            rho = 1000.0
            rho_e = Expression("rho", rho=rho, degree = 0)

            # loading, gravity in this case
            g = 9.81

            u_p_ = Function(U)
            u = TrialFunction(U)
            v = TestFunction(U)

            u_t, p_t = split(v)
            u_, p_ = split(u_p_)

            # kinematics with deformation gradient, right Green-Cauchy and Lagrangian strain tensor 
            I = Identity(3)
            F = variable(grad(u_) + I)
            C = variable(F.T * F)
            E = variable(0.5 * (C - I))

            # invariariants and Jacobian
            I_C = tr(C)
            II_C = (1.0 / 2.0) * (tr(C) ** 2 - tr(C * C))
            III_C = det(C)
            J = III_C ** (1.0 / 2.0)

            # stored strain energy density (incompressible Neo-Hooke model)
            psi = (Constant(mu) / 2.0) * (I_C - 2) - mu * ln(J) + p_ * ln(J) - (1.0 / (2.0 * Constant(lmbda))) * p_ ** 2

            # second Piola-Kirchoff conjugate with incremental Green-Lagrange
            S = 2.0 * diff(psi, C)
            D_E = derivative(E, u_p_, v)
            G = inner(S, D_E) * dx + inner((ln(J) - p_ / Constant(lmbda)), p_t) * dx - inner(Constant(g)*rho_e, -u_t[1])*dx
            J = derivative(G, u_p_, u)


            # define the ROI (Region Of Interest) corresponding to the bottom of the cube
            def bottom(x, on_boundary):
                return near(x[1], -1.0) and on_boundary


            # apply null Dirichlet boundary conditions on the bottom
            bcs = [DirichletBC(U.sub(0), Constant((0.0, 0.0, 0.0)), bottom)] 


            solver_parameters = {"newton_solver": {"linear_solver": "mumps", "maximum_iterations": 10}}
            try:
                solve(G == 0, u_p_, J=J, bcs=bcs, solver_parameters=solver_parameters)
            except RuntimeError:
                print("oups")

            # export the deformed mesh in the results directory
            ALE.move(mesh, u_p_.sub(0))
            with XDMFFile("results/forward_deformed_configuration_cube.xdmf") as xdmf:
                xdmf.write(mesh)


            ### INVERSE  ###

            # loading the deformed geometry
            mesh = Mesh()
            with XDMFFile("results/forward_deformed_configuration_cube.xdmf") as infile:
                infile.read(mesh)

            # mixed function space, 2nd order for displacement and 1st order for pressure
            V = VectorElement("CG", mesh.ufl_cell(), 2)
            P = FiniteElement("CG", mesh.ufl_cell(), 1)
            U = FunctionSpace(mesh, MixedElement([V, P]))

            # material properties
            mu = Constant(i)
            lmbda = Constant(j)
            rho = 1000.0
            rho_e = Expression("rho", rho=rho, degree = 0)

            # loading, gravity in this case
            g = 9.81

            u_p_ = Function(U)
            u = TrialFunction(U)
            v = TestFunction(U)

            u_t, p_t = split(v)
            u_, p_ = split(u_p_)

            # kinematics with deformation gradient, right Green-Cauchy and Lagrangian strain tensor 
            I = Identity(3)
            F = variable(inv(grad(u_) + I))
            C = variable(F.T * F)
            E = variable(0.5 * (C - I))

            # invariariants and Jacobian
            I_C = tr(C)
            II_C = (1.0 / 2.0) * (tr(C) ** 2 - tr(C * C))
            III_C = det(C)
            J = III_C ** (1.0 / 2.0)

            # stored strain energy density (incompressible Neo-Hooke model)
            psi = (Constant(mu) / 2.0) * (I_C - 2) - mu * ln(J) + p_ * ln(J) - (1.0 / (2.0 * Constant(lmbda))) * p_ ** 2

            # second Piola-Kirchoff conjugate with incremental Green-Lagrange
            S = 2.0*diff(psi, C)
            sigma = J**(-1.0)*F*S*F.T
            G = inner(sigma, sym(grad(u_t)))*dx + inner((ln(J) - p_/Constant(lmbda)), p_t)*dx - inner(Constant(g)*rho_e, -u_t[1])*dx
            J = derivative(G, u_p_, u)


            # define the ROI (Region Of Interest) corresponding to the top and bottom of the square
            def bottom(x, on_boundary):
                return near(x[1], -1.0) and on_boundary


            # apply null Dirichlet boundary conditions on the bottom
            bcs = [DirichletBC(U.sub(0), Constant((0.0, 0.0, 0.0)), bottom)] 

            solver_parameters = {"newton_solver": {"linear_solver": "mumps", "maximum_iterations": 10}}
            try:
                solve(G == 0, u_p_, J=J, bcs=bcs, solver_parameters=solver_parameters)
            except RuntimeError:
                print("oups")

            # export the deformed mesh in the results directory
            ALE.move(mesh, u_p_.sub(0))
            with XDMFFile("results/inverse_configuration_cube.xdmf") as xdmf:
                xdmf.write(mesh)


            ### ERROR  ###

            initial_mesh = Mesh()
            with XDMFFile("results/initial_configuration_cube.xdmf") as infile:
                infile.read(initial_mesh)

            inverse_mesh = Mesh()
            with XDMFFile("results/inverse_configuration_cube.xdmf") as infile:
                infile.read(inverse_mesh)

            initial_mesh_coordinates = initial_mesh.coordinates()
            inverse_mesh_coordinates = inverse_mesh.coordinates()

            mse = np.square(initial_mesh_coordinates - inverse_mesh_coordinates).mean()
            total_mse.append(mse)
            param_1.append(i)
            param_2.append(j)

    with open('results.npy', 'wb') as f:
        np.save(f, np.array([np.array(param_1), np.array(param_2), np.array(total_mse)]))

    print("export")


# Main body
if __name__ == '__main__':
    main()



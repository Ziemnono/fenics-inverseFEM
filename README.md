# fenics-inverseFEM
<h1>FEniCS code for inverse FEM hyper-elasticity</h1>
 
In this work, we developed a simple formulation and robust solution procedure using the FEniCS Project software for inverse deformation problems in soft-tissue biomechanics. Conversely to iterative algorithms, our method can solve with one simulation the rest-position without computing multiple solutions of the forward problem. Given a convergence threshold, our physics-based algorithm is about ten times faster and handles better large deformations than the iterative one.

Moreover, the approach also computes internal stresses in the organ with a one-shot approach, which does not require any additional direct deformation simulation from the rest configuration. The framework is implemented within an open-source pipeline enabling the seamless, fully parallelized, direct, and inverse deformation simulation of organs directly from segmented images. The design of the pipeline is flexible to user’s needs: for example, it allows the modification of the constitutive models by changing one single line of code.

# fenics-inverseFEM
<h2>FEniCS code for inverse FEM hyper-elasticity</h2>
 
In this work, we develop a framework for solving inverse deformation problems using the FEniCS Project finite-element software. We validate our approach with experimental imaging data acquired from a soft silicone beam under gravity. In contrast with inverse iterative algorithms that require multiple solutions of a standard elasticity problem, the proposed method can compute the undeformed configuration by solving only one modified elasticity problem. This modified problem has a complexity comparable to the standard one. The framework is implemented within an open-source pipeline enabling the direct and inverse deformation simulation directly from imaging data. We use the high-level unified form language (UFL) of the FEniCS Project to express the finite-element model in variational form and to automatically derive the consistent Jacobian. Consequently, the design of the pipeline is flexible: for example, it allows the modification of the constitutive models by changing a single line of code. We include a complete working example showing the inverse deformation of a beam deformed by gravity as supplementary material.


<h2>FEniCS</h2>
To download FEniCS: https://fenicsproject.org/download/
and know more about it: https://fenicsproject.org/

<h2>SOFA</h2>
We also used SOFA for comparing our FEniCS simulations.
To download SOFA: https://www.sofa-framework.org/download/ and know more about it: https://www.sofa-framework.org/.

We used the Multiplicative Jacobian Energy Decomposition (MJED) formulation, which is available in a specific SOFA plugin: https://github.com/sofa-framework/SofaMJEDFEM

<h2>Strcture of the repository</h2>
<h3>Submodule: FEgen</h3>
We used a custom code FEgen to generate tetrahedral volume meshes readbale by FEniCS (.xml) from surface meshes (.obj).
This tool can also read medical images NIFTI (.nii, .nii.gz) or INR (.inr, .inr.gz) and output VTK meshes (.vtu).

https://bitbucket.org/unilucompmech/fegen

<h3>example</h3>
The directory contains the different examples used in the paper. They have been implemented using the FEniCS Project software:

* Beam: hyperelastic cantilever beam deformed under gravity.
* Cube_compression: test case
* Cube_generalized_shear: appendix of the paper
* Cube_simple_shear: appendix of the paper
* meshes: containing a low resolution mesh of the beam used in the Beam repository. It also includes a python script for converting a .msh to .xdmf.
* Sagging_block: in section numerical results of the paper
 
<h3>sofa_scenes</h3>
The directory contains the sofa scenes used in the paper. You first need to compile SOFA with the MJED plugin:
* inData: contains different mesh refinement of the beam
* outData: results of each simulations
* scenes: the .py file used for creating and running the simulation

## Running the code
To simply run the codes:
```
python3 file_name.py
```

FEniCS also offers parallel computation:
```
mpirun -n N python3 file_name.py
```
where N is the number of core you want to use.

## Full paper
https://link.springer.com/article/10.1007/s00366-021-01597-z

## Information
### Authors 
- Arnaud Mazier: Department of Computational Science, Université du Luxembourg, Esch-sur-Alzette, Luxembourg
- Alexandre Bilger: Department of Computational Science, Université du Luxembourg, Esch-sur-Alzette, Luxembourg
- Antonio E. Forte: Harvard University,29 Oxford St,Cambridge MA 02138,USA
- Igor Peterlik: Institute of Computer Science, Masaryk University
- Jack S. Hale: Department of Computational Science, Université du Luxembourg, Esch-sur-Alzette, Luxembourg
- Stéphane P.A. Bordas: Department of Computational Science, Université du Luxembourg, Esch-sur-Alzette, Luxembourg

### Contact 
mazier.arnaud@gmail.com or use the [Discussion feature from Github](https://github.com/Ziemnono/fenics-inverseFEM/discussions)


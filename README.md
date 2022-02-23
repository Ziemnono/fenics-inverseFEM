# fenics-inverseFEM
<h2>FEniCS code for inverse FEM hyper-elasticity</h2>
 
In this work, we developed a simple formulation and robust solution procedure using the FEniCS Project software for inverse deformation problems in soft-tissue biomechanics. Conversely to iterative algorithms, our method can solve with one simulation the rest-position without computing multiple solutions of the forward problem. Given a convergence threshold, our physics-based algorithm is about ten times faster and handles better large deformations than the iterative one.

Moreover, the approach also computes internal stresses in the organ with a one-shot approach, which does not require any additional direct deformation simulation from the rest configuration. The framework is implemented within an open-source pipeline enabling the seamless, fully parallelized, direct, and inverse deformation simulation of organs directly from segmented images. The design of the pipeline is flexible to user’s needs: for example, it allows the modification of the constitutive models by changing one single line of code.


<h2>FEniCS</h2>
To download FEniCS: https://fenicsproject.org/download/
and know more about it: https://fenicsproject.org/

## Running the code
To simply run the code:
```
python3 file_name.py
```

FEniCS also offers parallel computation:
```
mpirun -n N python3 file_name.py
```
where N is the number of core you want to use.

<h2>Submodule: FEgen</h2>
We used a custom code FEgen to generate tetrahedral volume meshes readbale by FEniCS (.xml) from surface meshes (.obj).
This tool can also read medical images NIFTI (.nii, .nii.gz) or INR (.inr, .inr.gz) and output VTK meshes (.vtu).

https://bitbucket.org/unilucompmech/fegen

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
mazier.arnaud@gmail.com


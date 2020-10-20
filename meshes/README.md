To convert from gmsh file to xdmf file (widely used for FEniCS simulations) you can use: converter_from_msh_to_xdmf.py
Sometimes, we do obtain the following error:

  File "converter_from_msh_to_xdmf.py", line 15, in <module>
    for key in msh.cell_data_dict["gmsh:physical"].keys():
KeyError: 'gmsh:physical'

If you encounter this error you can use the command:
```
meshio-convert mesh_file.msh mesh_file.xdmf
```

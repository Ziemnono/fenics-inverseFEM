import meshio

filename = "Beam92k"

msh = meshio.read(filename + '.msh')

for cell in msh.cells:
    if cell.type == "triangle":
        triangle_cells = cell.data

    elif  cell.type == "tetra":
        tetra_cells = cell.data


for key in msh.cell_data_dict["gmsh:physical"].keys():

    if key == "triangle":
        triangle_data = msh.cell_data_dict["gmsh:physical"][key]

    elif key == "tetra":
        tetra_data = msh.cell_data_dict["gmsh:physical"][key]

tetra_mesh = meshio.Mesh(points=msh.points, cells={"tetra": tetra_cells})

meshio.write(filename + '.xdmf', tetra_mesh)


print("Conversion to xdmf done")
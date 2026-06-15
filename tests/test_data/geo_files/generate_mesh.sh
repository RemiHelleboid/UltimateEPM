# Script to generate the meshes from the geo files in the current directory.

list_2d_file=(square.geo disk.geo pn_2d.geo)
list_3d_file=(cube.geo sphere.geo)

for file in "${list_2d_file[@]}"; do
    gmsh -2 "$file" -format msh2 -o "../${file%.geo}.msh" -clmax 0.05
done

for file in "${list_3d_file[@]}"; do
    gmsh -3 "$file" -format msh2 -o "../${file%.geo}.msh" -clmax 0.075
done
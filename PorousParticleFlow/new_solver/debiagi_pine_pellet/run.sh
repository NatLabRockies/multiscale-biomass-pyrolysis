cp -r init 0
blockMesh
setFields
decomposePar
#mpirun -np 32 pyroFoam -parallel
pyroFoam

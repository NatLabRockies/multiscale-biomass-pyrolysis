#!/bin/sh
#cd ${0%/*} || exit 1    # Run from this directory

echo "copying init"
cp -r init 0
echo "running blockMesh"
blockMesh
echo "running snappyHexMesh"
snappyHexMesh -overwrite
echo "setting fields"
setFields
echo "setting porosity"
setPorosityField
decomposePar

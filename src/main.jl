using BioStructures
path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta/1ABO.rosetta_0001.pdb"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full/1ABO.pdb"
simulated = read(path, PDB)
real =read(truePath,PDB)

prop(x::AbstractVector{Bool}) = count(x)/length(x)

prop(resname.(simulated["A"]) .== resname.(real["A"]))


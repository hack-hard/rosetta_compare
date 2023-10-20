using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie

path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta/"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full/"
type=["1ABO","1BM2","1CKA","1G9O","1M61","1O4C","1R6J","2BYG"]

prop(x::AbstractVector{Bool}) = count(x)/length(x)

equality(x,y)=prop(resname.(x["A"]) .== resname.(y["A"]))
readdir()
function compare(type::String)::Tuple{Ref{ProteinStructure},Vector{ProteinStructure}}
    simulated = read.(path*"/".*filter(x ->(!isnothing∘match)(type*r".*rosetta.*pdb",x),readdir(path)), [PDB])
    real =read("$truePath/$type.pdb",PDB)
    (Ref(real),simulated)
end
type*r"\.*rosetta.*pdb"
scatter
σ(x) = std(x)/mean(x)

f(type::String) = equality.(compare(type)...)
f.(type)
mean(mean.(f.(type)))
std.(f.(type))
σ.(f.(type))
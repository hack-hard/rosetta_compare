using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie
using BioAlignments
using BioSymbols

path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta/"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full/"
type=["1ABO","1BM2","1CKA","1G9O","1M61","1O4C","1R6J","2BYG"]

BLOSUM62

prop(x::AbstractVector{Bool}) = count(x)/length(x)

equality(x,y)=prop(resname.(x["A"]) .== resname.(y["A"]))

function compare(type::String)::Tuple{Ref{ProteinStructure},Vector{ProteinStructure}}
    simulated = read.(path*"/".*filter(x ->(!isnothing∘match)(type*r".*rosetta.*pdb",x),readdir(path)), [PDB])
    real =read("$truePath/$type.pdb",PDB)
    (Ref(real),simulated)
end
Residue
AminoAcid
x = collect(compare(type[1])[1][]["A"])[1]
Base.convert(::Type{AminoAcid},r::Residue)= parse(AminoAcid,r.name)
σ(x) = std(x)/mean(x)
similarity(x::Residue,y::Residue) = BLOSUM62[x,y]
similarity(x::Chain,y::Chain) = mean(similarity.(x,y))
similarity(x::ProteinStructure,y::ProteinStructure) = similarity(x["A"],y["A"])
mean(similarity.(compare(type[1])...))
i(type::String) = equality.(compare(type)...)
s(type::String) = similarity.(compare(type)...)
i.(type)
mean(mean.(i.(type)))
mean(mean.(s.(type)))
std.(i.(type))
σ.(i.(type))
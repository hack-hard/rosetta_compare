using Revise
using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie
using BioAlignments
using BioSymbols
using CSV
using DataFrames
using BioSequences
using OffsetArrays

include("prot.jl")

path = "/home/hacquard/fullseqdesign_rosetta/allpositions"
rosetta_path = "$path/rosetta"
proteinMPNN_path = "$path/proteinMPNN"
reference_path = "/home/hacquard/fullseqdesign_rosetta/pdb.full"
protein = "/home/hacquard/rosetta_compare/data/protein.dat"
type = ["1ABO", "1CSK", "1CKA", "1R6J", "1G9O", "2BYG", "1BM2", "1O4C", "1M61"]
scores = ["score12", "talaris2013", "beta_nov16"]


p = CSV.read(protein, DataFrame; delim=" ", ignorerepeated=true)

p = hcat(p[!, 1:7], apply(l -> parse.(Int, split(l, ",")), p[!, 8:12]))




m = reshape(hcat(metric(equality), metric(similarity)), 1, :)
header = [:identity, :core_identity, :surface_identity, :similarity, :core_similarity, :surface_similarity]

res = compute.(m, scores) .|> Mesure
res = DataFrame(res, header)
res = hcat(DataFrame([scores], [:scores]), res)
CSV.write("temp.csv", res; transform=transform)
accessibility(f, p.corelist)(1)((j, j))
j = convert.(Chain, reference.(type))
transform(col, val) = val
transform(col, val::Mesure) = transform(col, val.mean) * " ± " * transform(col, val.std)
accessibility.(p.corelist, j)

# score12, talaris2013 proteinMPNN, beta_nov16
simulated(first(type), Val(:proteinMPNN))
dna(s::Model) = s |> only |> only |> collect .|> convert(AminoAcid)
function Base.read(io::IO, ::Type{Vector{Vector{AminoAcid}}})
    res = Vector{AminoAcid}[]
    while !eof(io)
        readline(io)
        push!(res, parse.(AminoAcid, collect(readline(io))))
    end
    res
end

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



p = CSV.read(protein, DataFrame; delim=" ", ignorerepeated=true)
p = hcat(p[!, 1:7], apply(l -> parse.(Int, split(l, ",")), p[!, 8:12]))


reference(first(type)) |> sequence

m = reshape(hcat(metric(equality), metric(similarity)), 1, :)
header = [:identity, :core_identity, :surface_identity, :similarity, :core_similarity, :surface_similarity]

res = compute.(m, Val(:rosetta),scores[1]) .|> Mesure
res = DataFrame(res, header)
res = hcat(DataFrame([scores], [:scores]), res)
CSV.write("temp.csv", res; transform=transform)
accessibility(f, p.corelist)(1)((j, j))
j = convert.(Chain, reference.(type))
transform(col, val) = val
transform(col, val::Mesure) = transform(col, val.mean) * " ± " * transform(col, val.std)
accessibility.(p.corelist, j)

# score12, talaris2013 proteinMPNN, beta_nov16

equality.(type |> first |> reference |> sequence |> parent |> Ref, simulated(first(type), Val(:proteinMPNN))[2:end])



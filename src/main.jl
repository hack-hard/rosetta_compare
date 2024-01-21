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





reference(first(type))
m = reshape(hcat(metric(equality), metric(similarity)), 1, :)
header = [:identity, :core_identity, :surface_identity, :similarity, :core_similarity, :surface_similarity]
simulated(1, Val(:proteinMPNN)) 
args = vcat(tuple.(Val(:rosetta),scores),(Val(:proteinMPNN),))


res = compute.(m, args) .|> Mesure
res = DataFrame(res, header)
column(values,header)= DataFrame([values],[header])
res = hcat(column(vcat(scores ,"proteinMPNN"), :scores), res)
CSV.write("temp.csv", res; transform=f)

f(col, val) = val
f(col, val::Mesure) = repr(val; context = :compact => true)
accessibility.(p.corelist, j)

# score12, talaris2013 proteinMPNN, beta_nov16

equality.(type |> first |> reference |> sequence |> parent |> Ref, simulated(first(type), Val(:proteinMPNN))[2:end])



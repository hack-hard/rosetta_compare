using Revise
using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie
using BioAlignments
using BioSymbols
using CSV
using DataFrames
using Printf
using BioSequences
path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full"
protein = "/home/hacquard/rosetta_compare/data/protein.dat"
type = ["1ABO", "1CSK", "1CKA", "1R6J", "1G9O", "2BYG","1BM2", "1O4C", "1M61"]
scores = ["score12", "talaris2013", "beta_nov16"]
Base.keys(s::ProteinStructure) = keys(s.models)
Base.keys(m::Model) = keys(m.chains)
Base.keys(c::Chain) = keys(c.residues)
Base.keys(r::Residue) = keys(r.atoms)

Base.convert(::Type{Chain}, m::Model) = only(m)
Base.convert(::Type{Model}, m::ProteinStructure) = only(m)
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)
Base.convert(::Type{Chain}, m::ProteinStructure) = convert(Chain, convert(Model, m))
prop(x::AbstractVector{Bool}) = count(x) / length(x)

equality(matrix::AbstractMatrix{Residue}) = equality(matrix[:, 1], matrix[:, 2])
equality(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = prop(equality.(x, y))
equality(x::Residue, y::Residue) = resname(x) == resname(y)

similarity(matrix::AbstractMatrix{Residue}) = similarity(matrix[1, :], matrix[2, :])
similarity(x::Residue, y::Residue) = BLOSUM62[x, y]
similarity(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = mean(similarity.(x, y))



simulated(type::String, score::String)::Vector{ProteinStructure} = (read.("$path/$score/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir("$path/$score/")), [PDB]))
reference(type::String)::ProteinStructure = read("$truePath/$type.pdb", PDB)

struct Partial <: Function
    f::Function
    args::Tuple
end
(f::Partial)(y...) = f.f(f.args..., y...)

partial(f, x...) = Partial(f, x)
pack(f::Function, x::Tuple) = f(x...)
pack(f) = partial(pack, f)


struct Mesure{T}
    mean::T
    std::T
end
Mesure(v::AbstractArray) = Mesure(mean(v), std(v))

function Base.show(io::IO, m::Mesure)
    show(io, m.mean)
    print(io, " ± ")
    show(io, m.std)
end

unpack(x) =
    if length(x) == 1
        only(x)
    else
        x
    end

multibroadcast(n::Integer, f, v...) =
    if n == 0
        f(v...)
    else
        unpack(multibroadcast.(n - 1, f, v...))
    end

σ(x) = std(x) / mean(x)

p = CSV.read(protein, DataFrame; delim=" ", ignorerepeated=true)
apply(f, d::DataFrame) = DataFrame(multibroadcast(2, f, eachcol(d)), names(d))
p = hcat(p[!, 1:7], apply(l -> parse.(Int, split(l, ",")), p[!, 8:12]))

(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)

compute(m::Function, score::String) = vcat([(multibroadcast(3, m(i), Ref.(reference(t)), simulated(t, score))) for (i, t) ∈ enumerate(type)]...)

function Base.filter(f, m::AbstractArray, dims::Integer)
    if dims > length(size(m))
        throw(ArgumentError("dims = $dims > $(length(size(m)))"))
    end
    stack(filter(f, eachslice(m; dims=dims)); dims=dims)
end

function Base.filter(f, v::Chain)
    l = filter(f, residues(v))
    Chain(v.id, collect(keys(l)), l, v.model)
end
accessibility(liste::AbstractVector{Int}, val::Chain) =
    filter(val) do (id, v)
        v.number ∈ liste
    end
accessibility(liste::AbstractVector{Int}, val::Chain...) = accessibility.(Ref(liste), val)
accessibility(f::Function, liste::AbstractVector{Int}, val::Chain...) = f(accessibility(liste, val...)...)
accessibility(f::Function, liste::AbstractVector{Int}) = partial(accessibility, f, liste)
accessibility(f::Function, liste::Vector{Vector{Int}}, i::Integer) = accessibility(f, liste[i])
accessibility(f::Function, liste::Vector{Vector{Int}}) = partial(accessibility, f, liste)
unif(f, x) = f
unif(f) = partial(unif, f)
metric(f::Function) = [unif(f), accessibility.(f, [p.corelist, p.surflist])...]
m = reshape(hcat(metric(equality), metric(similarity)), 1, :)
header = [:identity, :core_identity, :surface_identity, :similarity, :core_similarity, :surface_similarity]

f(x...) = x
res = compute.(m, scores) .|> Mesure
res = DataFrame(res,header)
res = hcat(DataFrame([scores],[:scores]),res)
CSV.write("temp.csv", res; transform=transform)
accessibility(f, p.corelist)(1)((j, j))
j = convert.(Chain, reference.(type))
transform(col,val) = val
transform(col,val::Mesure) = transform(col,val.mean)*" ± "* transform(col,val.std)
accessibility.(p.corelist, j)

# score12, talaris2013 proteinMPNN, beta_nov16

function dna(s::Chain)
function read_fa(io::IO)
    res = []
    while !eof(io)
        readline(io)
        append(res,parse(LongDNA,readline(io)))
    end
    res
end
        
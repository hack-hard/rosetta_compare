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
path = "/home/hacquard/fullseqdesign_rosetta/allpositions"
rosetta_path = "$path/rosetta"
proteinMPNN_path = "$path/proteinMPNN"
reference_path = "/home/hacquard/fullseqdesign_rosetta/pdb.full"
protein = "/home/hacquard/rosetta_compare/data/protein.dat"
type = ["1ABO", "1CSK", "1CKA", "1R6J", "1G9O", "2BYG","1BM2", "1O4C", "1M61"]
scores = ["score12", "talaris2013", "beta_nov16"]
Base.keys(s::ProteinStructure) = keys(s.models)
Base.keys(m::Model) = keys(m.chains)
Base.keys(c::Chain) = keys(c.residues)
Base.keys(r::Residue) = keys(r.atoms)
Base.convert(x) = partial(Base.convert,x)
Base.convert(::Type{Chain}, m::Model) = only(m)
Base.convert(::Type{Model}, m::ProteinStructure) = only(m)
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)
Base.convert(::Type{Chain}, m::ProteinStructure) = convert(Chain, convert(Model, m))


p = CSV.read(protein, DataFrame; delim=" ", ignorerepeated=true)
apply(f, d::DataFrame) = DataFrame(multibroadcast(2, f, eachcol(d)), names(d))
p = hcat(p[!, 1:7], apply(l -> parse.(Int, split(l, ",")), p[!, 8:12]))
p.firstdesignedpos
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
accessibility(liste::AbstractVector{Int}, val::AbstractVector) = val[liste]
accessibility(f::Function, liste::AbstractVector{Int}, val::AbstractVector) = f(accessibility(liste, val...)...)
accessibility(f::Function, liste::AbstractVector{Int}) = partial(accessibility, f, liste)
accessibility(f::Function, liste::Vector{Vector{Int}}, i::Integer) = accessibility(f, liste[i])
accessibility(f::Function, liste::Vector{Vector{Int}}) = partial(accessibility, f, liste)
unif(f, x) = f
unif(f) = partial(unif, f)
metric(f::Function) = [unif(f), accessibility.(f, [p.corelist, p.surflist])...]
m = reshape(hcat(metric(equality), metric(similarity)), 1, :)
header = [:identity, :core_identity, :surface_identity, :similarity, :core_similarity, :surface_similarity]

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
simulated(first(type),Val(:proteinMPNN))
dna(s::Model) = s|> only |> only |> collect .|> convert(AminoAcid)
function Base.read(io::IO,::Type{Vector{Vector{AminoAcid}}})
    res = Vector{AminoAcid}[]
    while !eof(io)
        readline(io)
        push!(res,parse.(AminoAcid,collect(readline(io))))
    end
    res
end

using Revise
using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie
using BioAlignments
using BioSymbols
using CSV
using DataFrames

path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta/"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full/"
protein = "/home/hacquard/rosetta_compare/data/protein.dat"
type = ["1ABO" "1CSK" "1BM2" "1CKA" "1G9O" "1M61" "1O4C" "1R6J" "2BYG"]

Base.keys(s::ProteinStructure) = keys(s.models)
Base.keys(m::Model) = keys(m.chains)
Base.keys(c::Chain) = keys(c.residues)
Base.keys(r::Residue) = keys(r.atoms)

Base.convert(::Type{Chain}, m::Model) = only(m)
Base.convert(::Type{Model}, m::ProteinStructure) = only(m)
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)
Base.convert(::Type{Chain},m::ProteinStructure) = convert(Chain,convert(Model,m))
prop(x::AbstractVector{Bool}) = count(x) / length(x)


equality(matrix::AbstractMatrix{Residue}) = equality(matrix[:, 1], matrix[:, 2])
equality(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = prop(equality.(x, y))
equality(x::Residue, y::Residue) = resname(x) == resname(y)

similarity(matrix::AbstractMatrix{Residue}) = similarity(matrix[1, :], matrix[2, :])
similarity(x::Residue, y::Residue) = BLOSUM62[x, y]
similarity(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = mean(similarity.(x, y))



simulated(type::String)::Vector{ProteinStructure} = (read.(path * "/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir(path)), [PDB]))[1:40]
reference(type::String)::ProteinStructure = read("$truePath/$type.pdb", PDB)

partial(f,x...)= (y...) -> f(x...,y...)
pack(f::Function, x::Tuple) = f(x...)
pack(f) = partial(pack,f)

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

p = CSV.read(protein,DataFrame; delim=" ", ignorerepeated=true)
apply(f,d::DataFrame) = DataFrame(multibroadcast(2,f , eachcol(d)), names(d))
p = hcat(p[!,1:7],apply(l -> parse.(Int,split(l,",")),p[!,8:12]))
p.hydrophcorelist
p.surflist
(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)

metric(m::Function) = [mean(multibroadcast(3, m(i), Ref.(reference(t)), simulated(t))) for (i,t) ∈ enumerate(type)]

function Base.filter(f, m::AbstractArray, dims::Integer)
    if dims > length(size(m))
        throw(ArgumentError("dims = $dims > $(length(size(m)))"))
    end
    stack(filter(f, eachslice(m; dims=dims)); dims=dims)
end

function Base.filter(f,v::Chain)
    l = filter(f,residues(v))
    Chain(v.id,collect(keys(l)),l,v.model)
end
accessibility(liste::AbstractVector{Int},val::Chain)  = filter(val) do (id,v)
    v.number ∈ liste
end
accessibility(liste::AbstractVector{Int},val::Chain ...)  = accessibility.(Ref(liste),val)

accessibility(f::Function,liste::AbstractVector{Int},val::Chain ...) = f(accessibility(liste,val...)...)
accessibility(f::Function) = partial(accessibility,f)
accessibility(f::Function, liste::AbstractVector{Int}) = partial(accessibility,f,liste)
accessibility(f::Function, liste::Vector{Vector{Int}}) = i -> accessibility(f,liste[i])


metric(accessibility(equality, p.surflist))
metric(accessibility(equality, p.corelist))

accessibility(pack(equality), p.surflist)(1)(i[1])
metric(similarity)
i = only(only(reference(type[4])))
j = only.(only.(simulated(type[4])))
accessibility.(equality,Ref(first(p.surflist)),Ref(i),j)
count.(x -> x  < 15 ,multibroadcast(2,accessibility,x))



y = only(only(first(simulated(first(type)))))


# score12, talaris2013 proteinMPNN, beta_nov16
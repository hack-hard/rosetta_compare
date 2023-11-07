using BioStructures
using Base.Filesystem
using Statistics
using CairoMakie
using BioAlignments
using BioSymbols

path = "/home/hacquard/fullseqdesign_rosetta/allpositions/rosetta/"
truePath = "/home/hacquard/fullseqdesign_rosetta/pdb.full/"
type = ["1ABO", "1BM2", "1CKA", "1G9O", "1M61", "1O4C", "1R6J", "2BYG"]

Base.keys(s::ProteinStructure) = keys(s.models)
Base.keys(m::Model) = keys(m.chains)
Base.keys(c::Chain) = keys(c.residues)
Base.keys(r::Residue) = keys(r.atoms)
Base.convert(::Type{Chain},m::Model) = only(m)
Base.convert(::Type{Model},m::ProteinStructure) = only(m)


prop(x::AbstractVector{Bool}) = count(x) / length(x)

equality(x::Chain, y::Chain) = prop(equality.(x,y))
equality(x::Residue, y::Residue) = resname(x) == resname(y)

function compare(type::String)::Tuple{Ref{ProteinStructure},Vector{ProteinStructure}}
    simulated = read.(path * "/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir(path)), [PDB])
    real = read("$truePath/$type.pdb", PDB)
    (Ref(real), simulated)
end


x = collect(only(compare(type[1])[1][]))[1]
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)
σ(x) = std(x) / mean(x)
similarity(x::Residue, y::Residue) = BLOSUM62[x, y]
similarity(x::Chain, y::Chain) = mean(similarity.(x, y))
similarity(x::ProteinStructure, y::ProteinStructure) = similarity(x["A"], y["A"])
mean(similarity.(compare(type[1])...))
i(type::String) = equality.(compare(type)...)
s(type::String) = similarity.(compare(type)...)
i.(type)
mean(mean.(i.(type)))
mean(mean.(s.(type)))
std.(i.(type))
σ.(i.(type))

(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)
accessibility(x::Atom) =x.temp_factor
accessibility(x::Residue) = mean(accessibility.(x))

unpack(x) = if length(x)==1
    only(x)
else
    x
end

multibroadcast(n::Integer, f,v...)= if n == 0
    f(v...)
else
    unpack(multibroadcast.(n-1,f,v...))
end


metric(f,ref::ProteinStructure,induit::ProteinStructure) =  multibroadcast(3, f,ref, induit)
metric.(Ref([equality, (x,_) -> accessibility(x)]),compare(first(type))... )[1]
metric.(Ref([equality]),compare(first(type))... )[1]

multibroadcast(2,x -> x+1,[[1],[2]])
ch = only(only(compare("1ABO")[1][]))

xs = rand(1:3, 1000)
ys = randn(1000)
violin(xs,ys)
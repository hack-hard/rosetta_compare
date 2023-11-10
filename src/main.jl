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

Base.convert(::Type{Chain}, m::Model) = only(m)
Base.convert(::Type{Model}, m::ProteinStructure) = only(m)
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)

prop(x::AbstractVector{Bool}) = count(x) / length(x)


equality(matrix::AbstractMatrix{Residue}) = equality(matrix[:, 1], matrix[:, 2])
equality(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = prop(equality.(x, y))
equality(x::Residue, y::Residue) = resname(x) == resname(y)

similarity(matrix::AbstractMatrix{Residue}) = similarity(matrix[1, :], matrix[2, :])
similarity(x::Residue, y::Residue) = BLOSUM62[x, y]
similarity(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = mean(similarity.(x, y))



simulated(type::String)::Vector{ProteinStructure} = (read.(path * "/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir(path)), [PDB]))[1:40]
reference(type::String)::ProteinStructure = read("$truePath/$type.pdb", PDB)



σ(x) = std(x) / mean(x)



(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)
accessibility(x::Atom) = x.temp_factor
accessibility(x::Residue) = mean(accessibility.(x))

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


metric(m) = mean.(multibroadcast(4, m, Ref.(reference.(type)), simulated.(type)))

function Base.filter(f, m::AbstractArray, dims::Integer)
    if dims > length(size(m))
        throw(ArgumentError("dims = $dims > $(length(size(m)))"))
    end
    stack(filter(f, eachslice(m; dims=dims)); dims=dims)
end

filter_accessibility(predicate::Function) = (ref::Chain, val::Chain) -> filter(hcat(collect(ref), collect(val)), 1) do r
    ref = r[1]
    println(typeof(ref))
    predicate(accessibility(ref))
end

core(metric) = metric ∘ filter_accessibility(x -> x < 15)
surface(metric) = metric ∘ filter_accessibility(x -> x > 0.3)

metric(core(equality))

metric(similarity)

x = only.(only.(reference.(type)))
count.(x -> x  < 15 ,multibroadcast(2,accessibility,x))

density(accessibility.(x[7]))

density(accessibility.(x))
y = only(only(first(simulated(first(type)))))


multibroadcast(4, equality ∘ filter_accessibility(x -> x < .15),Ref.(reference.(type)), simulated.(type)) .|> mean
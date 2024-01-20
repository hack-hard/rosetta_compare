prop(x::AbstractVector{Bool}) = count(x) / length(x)

equality(matrix::AbstractMatrix{Residue}) = equality(matrix[:, 1], matrix[:, 2])
equality(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = prop(equality.(x, y))
equality(x::Residue, y::Residue) = resname(x) == resname(y)

similarity(matrix::AbstractMatrix{Residue}) = similarity(matrix[1, :], matrix[2, :])
similarity(x::Residue, y::Residue) = BLOSUM62[x, y]
similarity(x::Union{Chain,AbstractVector{Residue}}, y::Union{Chain,AbstractVector{Residue}}) = mean(similarity.(x, y))


simulated(type::String,model::Val{:rosetta}, score::String)::Vector{Vector{AminoAcid}} = read.("$path/$score/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir("$rosetta_path/$score/")), [PDB]) .|> dna
simulated(type::String,model::Val{:proteinMPNN}) = read("$proteinMPNN_path/seqs/$type.fa" , Vector{Vector{AminoAcid}})
reference(type::String)::ProteinStructure = read("$reference_path/$type.pdb", PDB)

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


Base.keys(s::ProteinStructure) = keys(s.models)
Base.keys(m::Model) = keys(m.chains)
Base.keys(c::Chain) = keys(c.residues)
Base.keys(r::Residue) = keys(r.atoms)
Base.convert(x) = partial(Base.convert,x)
Base.convert(::Type{Chain}, m::Model) = only(m)
Base.convert(::Type{Model}, m::ProteinStructure) = only(m)
Base.convert(::Type{AminoAcid}, r::Residue) = parse(AminoAcid, r.name)
Base.convert(::Type{Chain}, m::ProteinStructure) = convert(Chain, convert(Model, m))


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

compute(m::Function, score::String) = vcat([(multibroadcast(3, m(i), Ref.(reference(t)), simulated(t, score))) for (i, t) ∈ enumerate(type)]...)


unif(f) = partial(unif, f)
metric(f::Function) = [unif(f), accessibility.(f, [p.corelist, p.surflist])...]

apply(f, d::DataFrame) = DataFrame(multibroadcast(2, f, eachcol(d)), names(d))
(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)

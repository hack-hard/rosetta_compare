path = "/home/hacquard/fullseqdesign_rosetta/allpositions"
rosetta_path = "$path/rosetta"
proteinMPNN_path = "$path/proteinMPNN"
reference_path = "/home/hacquard/fullseqdesign_rosetta/pdb.full"
protein = "/home/hacquard/rosetta_compare/data/protein.dat"
type = ["1ABO", "1CSK", "1CKA", "1R6J", "1G9O", "2BYG", "1BM2", "1O4C", "1M61"]
scores = ["score12", "talaris2013", "beta_nov16"]


prop(x::AbstractVector{Bool}) = count(x) / length(x)

equality(x::AbstractVector{AminoAcid}, y::AbstractVector{AminoAcid}) = prop(equality.(x, y))
equality(x::AminoAcid, y::AminoAcid) = x == y

similarity(x::AminoAcid, y::AminoAcid) = BLOSUM62[x, y]
similarity(x::AbstractVector{AminoAcid}, y::AbstractVector{AminoAcid}) = mean(similarity.(x, y))


simulated(type::String,model::Val{:rosetta}, score::String)::Vector{Vector{AminoAcid}} = read.("$rosetta_path/$score/" .* filter(x -> (!isnothing ∘ match)(type * r".*rosetta.*pdb", x), readdir("$rosetta_path/$score/")), [PDB]) .|> dna
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
Base.parse(x::Type) = partial(parse,x)

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
accessibility(f::Function, liste::AbstractVector{Int}, val::AbstractVector...) = f(accessibility.(Ref(liste), val)...)
accessibility(f::Function, liste::AbstractVector{Int}) = partial(accessibility, f, liste)
accessibility(f::Function, liste::Vector{Vector{Int}}, i::Integer) = accessibility(f, liste[i])
accessibility(f::Function, liste::Vector{Vector{Int}}) = partial(accessibility, f, liste)

sequence(s::ProteinStructure) = s |> only |> only |> sequence
sequence(s::Chain) = OffsetVector(s |> collect .|> convert(AminoAcid),first_pos(s)-1)

first_pos(s::Chain) =  s |> keys .|> parse(Int) |> minimum


unif(f, _) = f

compute(m::Function, generators...) = mapreduce(vcat,enumerate(type)) do (i,t)
    m(i).(reference(t) |> sequence |> Ref, simulated(t, generators...) .|> sequence)
end


unif(f) = partial(unif, f)
metric(f::Function) = [unif(f), accessibility.(f, [p.corelist, p.surflist])...]

apply(f, d::DataFrame) = DataFrame(multibroadcast(2, f, eachcol(d)), names(d))
(v::AbstractVector{Function})(x...) = map(f -> f(x...), v)

function Base.read(io::IO, ::Type{Vector{Vector{AminoAcid}}})
    res = Vector{AminoAcid}[]
    while !eof(io)
        readline(io)
        push!(res, parse.(AminoAcid, collect(readline(io))))
    end
    res
end
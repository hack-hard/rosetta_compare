using DataFrames
using CSV
trueData = "/home/hacquard/fullseqdesign_rosetta/pdb.full"
function loadpdb(file::String)::DataFrame
    DataFrame(CSV.File(file; header = [:atom, :id, :element,:radical,:chain,:radical_id,:x,:y,:z,:occupancy,:temperature]))
end
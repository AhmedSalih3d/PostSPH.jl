__precompile__()

module PostSPH

using Printf #To construct string message easily
using StaticArrays

#Add to project.toml manually: https://discourse.julialang.org/t/update-project-toml-manually/32477

include("ReadXML.jl")
include("SaveVTK.jl")

export
    Cat,
    readBi4Array,
    readBi4Body,
    readBi4Particles,
    readBi4Time,
    MkArray

##Hardcoded enum - Cat is "category"
@enum Cat begin
            Points
            Idp
            Vel
            Rhop
          end

##Read from binary files directly

## Make two dicts,   one for the enum -> string and one for the data type - Bi4 specific
const varNames     = (:key, :offset)

const IdpSearch    = transcode(UInt8,"ARRAY\x03\0\0\0Idp")[:]
const IdpOffset    = 11
const IdpKey       = NamedTuple{varNames}([IdpSearch,IdpOffset])

const PosSearch    = transcode(UInt8,"ARRAY\x03\0\0\0Pos")[:]
const PosOffset    = 11
const PosKey       = NamedTuple{varNames}([PosSearch,PosOffset])

const VelSearch    = transcode(UInt8,"ARRAY\x03\0\0\0Vel")[:]
const VelOffset    = 11
const VelKey       = NamedTuple{varNames}([VelSearch,VelOffset])

const RhopSearch   = transcode(UInt8,"ARRAY\x04\0\0\0Rhop")[:]
const RhopOffset   = 12
const RhopKey      = NamedTuple{varNames}([RhopSearch,RhopOffset])

const searchKeyBi4    = Dict{Cat,NamedTuple}(Idp => IdpKey, Points => PosKey, Vel => VelKey, Rhop => RhopKey)
const catTypeBi4      = Dict{Cat,DataType}(Idp => Int32, Points => Float32, Vel => Float32, Rhop => Float32)
const catArrayBi4     = Dict{Cat,Int64}(Idp => 1, Points => 2, Vel => 2, Rhop => 1)
const catColBi4       = Dict{Cat,Int64}(Idp => 1, Points => 3, Vel => 3, Rhop => 1)

##Lists files in directory and only returns applicable files, ie. "Part_XXXX.bi4"
# first_file bug
function _dirFiles(first_file::Bool=false)
    if first_file == false
        files = readdir()
    else
        files = readdir()[1]
    end
    #Operation on dirFiles instantly
    filter!(x->occursin(r"Part_\d{4}.bi4",x),files)
    return files
end

#Command for plotting
#println(Char.(key))
#scatter(getindex.(a[1], 1), getindex.(a[1], 3))
#readBi4Array(typ::Cat,Bi4Files::Array{String,1}=_dirFiles()) = readBi4Array(typ, false, Bi4Files)
#readBi4Array(typ::Cat,Bi4Files::Array{String,1}) = readBi4Array(typ, false, Bi4Files)
readBi4Array(typ::Cat,Bi4Files::String) = readBi4Array(typ, [Bi4Files])
function readBi4Array(typ::Cat,Bi4Files::Array{String,1}=_dirFiles())

    nBi4     = size(Bi4Files)[1]

    key    = searchKeyBi4[typ].key
    offset = searchKeyBi4[typ].offset
    ncol   = catColBi4[typ]
    T      = catTypeBi4[typ]

    # THIS BREAKS SAVE VTK
    j = Vector{Array{T,1}}(undef,nBi4)
    Threads.@threads for i = 1:nBi4
        j_tmp,~ = _readBi4(Bi4Files[i],key,offset,T,ncol)
        j[i] = j_tmp
    end

    return j
end


function _readBi4(file::String,key,offset,T,ncol)

    # Import a full bi4 file as Array{UInt8,1}
    rf = _rf(file)

    # Start position is found by search the file for the key and finding
    # first occurence, then adding offset
    startPos = Base._searchindex(rf, key, 1) + offset #1 byte offset

    #+1 due to hexeditor/julia?. "n id start" and "n id end"
    nid_s = startPos + 1  + sizeof(Int64)
    nid_e = nid_s    - 1  + sizeof(Int32)

    # Every "Array" in bi4 mentions the number of particles, read it in
    n = reinterpret(Int32,rf[nid_s:nid_e])[1]

    #data id start
    # Multiply with 4 here since UInt8 size, times number of particles, times
    # times number of columns gives the correct indices in the rf array for
    # Float32, Int32 etc.
    did_s = nid_e + 1 + sizeof(Int32)
    did_e = did_s - 1  + 4*n*ncol

    # Reinterpret the data as the specified data type, extract the relevant
    # snip of Array{UInt8,1} in "rf"
    data  = reinterpret(T,rf[did_s:did_e])

    return data,n
end

# Returns Uint8 array of file
function _rf(file::String)
    ft = open(file,read=true)
    rf = read(ft)
    close(ft)
    return rf
end
##StaticArrays is hard to use here since it is needed to offset with "Int32", between
# all searches unlike "readBi4Array"
# Useless since these values inside do not change
function readBi4Particles(Bi4Files::Array{String,1}=_dirFiles())

    nBi4     = size(Bi4Files)[1]

    j  = Vector{Array{Int32,1}}(undef,nBi4)

    ParticleString = ["CaseNp","CaseNfixed","CaseNmoving","CaseNfloat","CaseNfluid"]
    for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        j_inner = zeros(Int32,5)
        for k in eachindex(j_inner)
            readuntil(ft,ParticleString[k])
            read(ft,Int32)
            j_inner[k] = read(ft,Int32)
        end
        j[i] = j_inner
        close(ft)
    end
    return j
end

#Npok is the current number of actual particles in the bi4 file
function readBi4Npok(Bi4Files::Array{String,1}=_dirFiles())

    nBi4     = size(Bi4Files)[1]

    j  = Vector{Array{Int32,1}}(undef,nBi4)

    ParticleString = "Npok"
    j              = zeros(Int32,(nBi4,))
    for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        readuntil(ft,ParticleString)
        read(ft,Int32)
        j[i] = read(ft,Int32)
        close(ft)
    end
    return j
end

## Function to read the time at current simulation step, as in "XXXX.out"
function readBi4Time(Bi4Files::Array{String,1}=_dirFiles())
    nBi4     = size(Bi4Files)[1]

    T  = Float64
    j  = Vector{T}(undef,nBi4)

    ParticleString = ["TimeStep"]
    Threads.@threads for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        readuntil(ft,ParticleString[1])
        read(ft,Int32)
        j[i] = read(ft,T)
        close(ft)
    end
    return j
end

#Function to only find specific Idps
#for MovingSquare example "Bodies[2][1]""
#In the for loop the first index is the index of the relevant "typ" array,
#the second index is the sorting of this "typ" array corresponding to the "idp"
#array and "start:move" are the number of particles defined by the "Body", where
#0 and 1 indexing from C++ and Julia differences has been taken into account.
function readBi4Body(Body,typ)
    start = getfield(Body, :beg)+1    #First idp

    idp_vec  = readBi4Array(Idp)
    val_vec  = readBi4Array(typ)

    nBi4     = length(idp_vec)

    k = []
    move = []

    if Body.bool == true
        move  = readBi4Npok( )
        k     = collect(1:nBi4)
    else
        move  = start+getfield(Body, :count)-1  #Number of particles from first idp
        k     = ones(Int,nBi4)
    end


    j = similar(val_vec)

    if typ == Idp || typ == Rhop
        for i = 1:length(j)
            j[i] = val_vec[i][sortperm(idp_vec[i])][start:move[k[i]]]
        end
    else
        for i = 1:length(j)
            id = sortperm(idp_vec[i])
            j[i] = val_vec[i][:,id][:,start:move[k[i]]]
        end
    end

    return j
end

end #PostSPH

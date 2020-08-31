__precompile__()

module PostSPH2

using Printf #To construct string message easily
using StaticArrays

include("SaveVTK.jl")

export
    readVtkArray,
    Cat,
    readVtkParticles,
    ForceVtk,
    FloatingVtkTranslation,
    MassVtk,
    readVtkNames,
    readVtkVariables,
    readBi4Array,
    readBi4Body,
    MkArray

##Hardcoded enum - Cat is "category"
        @enum Cat begin
            Points
            Idp
            Vel
            Rhop
            Mass
            Press
            Vol
            Ace
            Vor
            Typ #Typ, since "Type" is illegal to assign
            Mk
        end

##Read from binary files directly

## Make two dicts,   one for the enum -> string and one for the data type - Bi4 specific
const varNames  = (:key, :offset)

const IdpSearch    = [0x41;0x52;0x52;0x41;0x59;0x03;0x00;0x00;0x00;0x49;0x64;0x70]
const IdpOffset    = 11
const IdpKey       = NamedTuple{varNames}([IdpSearch,IdpOffset])

const PosSearch    = [0x41;0x52;0x52;0x41;0x59;0x03;0x00;0x00;0x00;0x50;0x6f;0x73]
const PosOffset    = 11
const PosKey       = NamedTuple{varNames}([PosSearch,PosOffset])

const VelSearch    = [0x41;0x52;0x52;0x41;0x59;0x03;0x00;0x00;0x00;0x56;0x65;0x6c]
const VelOffset    = 11
const VelKey       = NamedTuple{varNames}([VelSearch,VelOffset])

const RhopSearch   = [0x41;0x52;0x52;0x41;0x59;0x04;0x00;0x00;0x00;0x52;0x68;0x6f;0x70]
const RhopOffset   = 12
const RhopKey      = NamedTuple{varNames}([RhopSearch,RhopOffset])

const searchKeyBi4    = Dict{Cat,NamedTuple}(Idp => IdpKey, Points => PosKey, Vel => VelKey, Rhop => RhopKey)
const catTypeBi4 = Dict{Cat,DataType}(Idp => Int32, Points => Float32, Vel => Float32, Rhop => Float32)
const catArrayBi4  = Dict{Cat,Int64}(Idp => 1, Points => 2, Vel => 2, Rhop => 1)
const catColBi4  = Dict{Cat,Int64}(Idp => 1, Points => 3, Vel => 3, Rhop => 1)

##Lists files in directory and only returns applicable files, ie. "Part_XXXX.bi4"
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

##transferData is smart since we change the value of an array in place and
# skip a lot of unnecessary allocation steps. It has been tuned to accept
# matrices and vectors, where the type is automatically deduced using "eltype".
 _readel(ft::IOStream, typ) = read(ft, typ)
@inline _readel(ft::IOStream, typ::Type{<:SVector{3,T}}) where T = SVector(_readel(ft, T), _readel(ft, T), _readel(ft, T))

function _transferDataBi4(ft::IOStream, arrayVal::AbstractVector)
    typ = eltype(arrayVal)
    @inbounds for i in eachindex(arrayVal)
        arrayVal[i] = _readel(ft, typ)
    end
end


##If nCol is not 1, then type is static array
_typeMaker(typ, nCol) = nCol == 1 ? typ : SVector{nCol, typ}

#Command for plotting
#println(Char.(key))
#scatter(getindex.(a[1], 1), getindex.(a[1], 3))
#readBi4Array(typ::Cat,Bi4Files::Array{String,1}=_dirFiles()) = readBi4Array(typ, false, Bi4Files)
readBi4Array(typ::Cat,Bi4Files::Array{String,1}) = readBi4Array(typ, false, Bi4Files)
readBi4Array(typ::Cat,SeekNull::Bool,Bi4Files::String) = readBi4Array(typ, false, [Bi4Files])
readBi4Array(typ::Cat,Bi4Files::String) = readBi4Array(typ, false, [Bi4Files])
function readBi4Array(typ::Cat,SeekNull::Bool=false,Bi4Files::Array{String,1}=_dirFiles())

    #if true start from top, else do not, use startPos from _Bi4Pos
    nBi4     = size(Bi4Files)[1]

    T  = _typeMaker(catTypeBi4[typ], catColBi4[typ])
    j  = fill(Array{T,1}(), nBi4) #Less allocs than Vector{Array{T}}(undef,nBi4)

    key    = searchKeyBi4[typ].key
    offset = searchKeyBi4[typ].offset

    Threads.@threads for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        rf = read(ft)
        startPos = Base._searchindex(rf, key, 1) + offset #1 byte offset
        seek(ft,startPos)

            read(ft,Int64)
        n = read(ft,Int32)
            read(ft,Int32)

            j[i] = zeros(T,n)
            _transferDataBi4(ft,j[i])


        close(ft)
    end
    return j
end

##StaticArrays is hard to use here since it is needed to offset with "Int32", between
# all searches unlike "readBi4Array"
function readBi4Particles(Bi4Files::Array{String,1}=_dirFiles())

    nBi4     = size(Bi4Files)[1]

    j  = Vector{Array{Int32,1}}(undef,nBi4)

    ParticleString = ["CaseNp","CaseNfixed","CaseNmoving","CaseNfloat","CaseNfluid"]
    Threads.@threads for i = 1:nBi4
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
readBi4Body(Body,typ,Bi4Files::String) =  readBi4Body(Body, typ, [Bi4Files])
function readBi4Body(Body,typ,Bi4Files::Array{String,1}=_dirFiles())
    start = getfield(Body, :beg)+1    #First idp
    move  = start+getfield(Body, :count)-1  #Number of particles from first idp
    idp_vec  = readBi4Array(Idp,false,Bi4Files)
    val_vec  = readBi4Array(typ,false,Bi4Files)

    j = similar(val_vec)

    for i = 1:length(j)
        j[i] = val_vec[i][sortperm(idp_vec[i])][start:move]
    end

    return j
end

#Functions



end #PostSPH2

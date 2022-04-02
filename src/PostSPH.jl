__precompile__()

module PostSPH
using OrderedCollections

#Add to project.toml manually: https://discourse.julialang.org/t/update-project-toml-manually/32477
include("SaveVTK.jl")

export Cat,
    readBi4Array, readBi4_NumberOfParticles, readBi4_CurrentTotalParticles, readBi4_Time, readBi4_Head

##Hardcoded enum - Cat is "category"
"""
    Cat are hardcoded enums to specify the potential arrays extractable from bi4 files.
"""
@enum Cat begin
    Points
    Idp
    Vel
    Rhop
end

##Read from binary files directly

## Make two dicts,   one for the enum -> string and one for the data type - Bi4 specific
const varNames = (:key, :offset)

const IdpSearch = transcode(UInt8, "ARRAY\x03\0\0\0Idp")[:]
const IdpOffset = 11
const IdpKey = NamedTuple{varNames}([IdpSearch, IdpOffset])

const PosSearch = transcode(UInt8, "ARRAY\x03\0\0\0Pos")[:]
const PosOffset = 11
const PosKey = NamedTuple{varNames}([PosSearch, PosOffset])

const VelSearch = transcode(UInt8, "ARRAY\x03\0\0\0Vel")[:]
const VelOffset = 11
const VelKey = NamedTuple{varNames}([VelSearch, VelOffset])

const RhopSearch = transcode(UInt8, "ARRAY\x04\0\0\0Rhop")[:]
const RhopOffset = 12
const RhopKey = NamedTuple{varNames}([RhopSearch, RhopOffset])

const searchKeyBi4 = Dict{Cat,NamedTuple}(Idp => IdpKey, Points => PosKey, Vel => VelKey, Rhop => RhopKey)
const catTypeBi4   = Dict{Cat,DataType}(Idp => Int32, Points => Float32, Vel => Float32, Rhop => Float32)
const catColBi4    = Dict{Cat,Int64}(Idp => 1, Points => 3, Vel => 3, Rhop => 1)

##Lists files in directory and only returns applicable files, ie. "Part_XXXX.bi4"
function _dirFiles(rgxPat::Regex=Regex("Part_\\d{4}.bi4"))
    files = readdir()
    #Operation on dirFiles instantly
    filter!(x -> occursin(rgxPat, x), files)
    return files
end


readBi4Array(typ::Cat, Bi4Files::String) = readBi4Array(typ, [Bi4Files])
function readBi4Array(typ::Cat, Bi4Files::Vector{String} = _dirFiles())

    key    = searchKeyBi4[typ].key
    offset = searchKeyBi4[typ].offset
    T      = catTypeBi4[typ]
    ncol   = catColBi4[typ]
   
    j  = similar(Vector{Vector{T}},axes(Bi4Files))

    Base.Threads.@threads for i in eachindex(j)
        println("$typ | Iteration: " * lpad(string(i), 4, "0") * "|Reading: " * Bi4Files[i])
        j[i], ~ = _readBi4(Bi4Files[i], key, offset, T, ncol)
    end

    return j
end


function _readBi4(file::String, key::Vector{UInt8}, offset::Int, T::DataType, ncol::Int)

    # Import a full bi4 file as Array{UInt8,1}
    ft = open(file, read = true)
    rf = read(ft)
    close(ft)

    # Start position is found by search the file for the key and finding
    # first occurence, then adding offset
    startPos = Base._searchindex(rf, key, 1) + offset #1 byte offset

    #+1 due to hexeditor/julia?. "n id start" and "n id end"
    nid_s = startPos + 1 + sizeof(Int64)
    nid_e = nid_s - 1 + sizeof(Int32)

    # Every "Array" in bi4 mentions the number of particles, read it in
    n = reinterpret(Int32, rf[nid_s:nid_e])[1]

    #data id start
    # Multiply with 4 here since UInt8 size, times number of particles, times
    # times number of columns gives the correct indices in the rf array for
    # Float32, Int32 etc.
    did_s = nid_e + 1 + sizeof(Int32)
    did_e = did_s - 1 + 4 * n * ncol

    # Reinterpret the data as the specified data type, extract the relevant
    # snip of Array{UInt8,1} in "rf"
    data = reinterpret(T, rf[did_s:did_e])

    return data, n
end


function readBi4_NumberOfParticles(Bi4Files::Vector{String} = _dirFiles())

    j  = similar(Vector{Vector{Int32}},axes(Bi4Files))

    ParticleString = ["CaseNp", "CaseNfixed", "CaseNmoving", "CaseNfloat", "CaseNfluid"]
    for i in eachindex(j)
        ft = open(Bi4Files[i], read = true)
        j_inner = zeros(Int32, length(ParticleString))
        for k in eachindex(j_inner)
            readuntil(ft, ParticleString[k])
            read(ft, Int32)
            j_inner[k] = read(ft, Int32)
        end
        j[i] = j_inner
        close(ft)
    end
    return ParticleString, j
end

#Npok is the current number of actual particles in the bi4 file
function readBi4_CurrentTotalParticles(Bi4Files::Vector{String} = _dirFiles())

    nBi4 = size(Bi4Files)[1]

    ParticleString = "Npok"
    j = zeros(Int32, (nBi4,))
    for i = 1:nBi4
        ft = open(Bi4Files[i], read = true)
        readuntil(ft, ParticleString)
        read(ft, Int32)
        j[i] = read(ft, Int32)
        close(ft)
    end
    return j
end

## Function to read the time at current simulation step, as in "XXXX.out"
function readBi4_Time(Bi4Files::Vector{String} = _dirFiles())
    nBi4 = size(Bi4Files)[1]

    T = Float64
    j = Vector{T}(undef, nBi4)

    ParticleString = ["TimeStep"]
    Threads.@threads for i = 1:nBi4
        ft = open(Bi4Files[i], read = true)
        readuntil(ft, ParticleString[1])
        read(ft, Int32)
        j[i] = read(ft, T)
        close(ft)
    end
    return j
end

function readBi4_Head()
    Bi4Head = _dirFiles(Regex("Part_Head"))

    file    = Bi4Head[1]

    # Import a full bi4 file as Array{UInt8,1}
    ft = open(file, read = true)
    rf = read(ft)
    close(ft)

    # Start position is found by search the file for the key and finding
    # first occurence, then adding offset
    key      = codeunits("ITEM")
    offset   = -1               #To start at "I" of "ITEM"

    # Bad code
    loc_val   = 1;
    Item_Locs = Vector{Int64}()
    while(true)
        #+2 to jump up to the next search ITEM!
        loc_val = Base._searchindex(rf, key, loc_val+2) + offset #1 byte offset
        if loc_val == -1
            break
        else
            push!(Item_Locs,loc_val)
        end
    end
    popfirst!(Item_Locs)          #Remove 1 info Item
    popfirst!(Item_Locs)          #Remove 2 info Item
    push!(Item_Locs,length(rf)+1) #Add last range

    Item_Ranges = Vector{UnitRange{Int64}}()
    for i = 1:length(Item_Locs)-1
        push!(Item_Ranges,Item_Locs[i]:(Item_Locs[i+1]-1))
    end

    
    function searchValue(str2Search::Vector{UInt8},strNeedle::String,seekCounter::Int,OutputType::DataType)
        loc_a = Base._searchindex(str2Search,codeunits(strNeedle),seekCounter)+ncodeunits(strNeedle)+sizeof(Int32)
        loc_b = loc_a+(sizeof(Int32)-1)
        range_ab = loc_a:loc_b

        valR      = reinterpret(OutputType,str2Search[range_ab])[1]

        return valR
    end

    function searchType(str2Search::Vector{UInt8})
        possibleTypes = ["Fixed";"Floating";"Fluid"]

        for pT in possibleTypes
            checkVal = Base._searchindex(str2Search,codeunits(pT),1)

            if sign(checkVal) == 0
                continue
            else
                return pT
            end
        end
    end

    Bi4_IdCount = 0
    dct = Vector{OrderedDict}(undef,length(Item_Ranges))
    for (ival,valRange) in enumerate(Item_Ranges)
        rf_ = rf[valRange]
        Count = searchValue(rf_,"Count",1,Int32)

        MkType = searchValue(rf_,"MkType",1,Int32)

        # Have to skip "MkBlocks" syntax..
        Mk = searchValue(rf_,"Mk",20,Int32)

        ActualType_ = searchType(rf_)

        IdRange = (Bi4_IdCount+1):(Bi4_IdCount+Count-1) #+1 -1 -> Julia Indexing
        Bi4_IdCount += Count;

        dct[ival] = OrderedDict("Type"=>ActualType_,"MkType"=>MkType,"Mk"=>Mk, "Count"=>Count, "IdRangeJulia"=>IdRange,"IdRangeBi4"=>(IdRange).-1)
    end

    return dct

end

end #PostSPH

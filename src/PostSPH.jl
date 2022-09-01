__precompile__()

module PostSPH
using OrderedCollections
using Printf

#Add to project.toml manually: https://discourse.julialang.org/t/update-project-toml-manually/32477
include("./SaveVTK.jl")

export Bi4DataType,
    readBi4Array, readBi4_NumberOfParticles, readBi4_CurrentTotalParticles, readBi4_Time, readBi4_Head

"""
    Bi4DataType are hardcoded enums to specify the potential arrays extractable from bi4 files.
"""
@enum Bi4DataType begin
    Points
    Points_F64
    Idp
    Vel
    Rhop
end

##Read from binary files directly

## Final Dict
const varNames = NamedTuple{(:key,:offset,:type,:ncol), Tuple{Vector{UInt8},Int,DataType,Int}}
_Bi4DataTypeDict             = Dict{Bi4DataType,varNames}()
_Bi4DataTypeDict[Points]     = (key = transcode(UInt8, "ARRAY\x03\0\0\0Pos")[:]  , offset = 11, type = Float32   , ncol   = 3 )
_Bi4DataTypeDict[Points_F64] = (key = transcode(UInt8, "ARRAY\x04\0\0\0Posd")[:] , offset = 12, type = Float64   , ncol   = 3 )
_Bi4DataTypeDict[Idp]        = (key = transcode(UInt8, "ARRAY\x03\0\0\0Idp")[:]  , offset = 11, type = Int32     , ncol   = 1 )
_Bi4DataTypeDict[Vel]        = (key = transcode(UInt8, "ARRAY\x03\0\0\0Vel")[:]  , offset = 11, type = Float32   , ncol   = 3 )
_Bi4DataTypeDict[Rhop]       = (key = transcode(UInt8, "ARRAY\x04\0\0\0Rhop")[:] , offset = 12, type = Float32   , ncol   = 1 )
const Bi4DataTypeDict        = deepcopy(_Bi4DataTypeDict)


##Lists files in directory and only returns applicable files, ie. "Part_XXXX.bi4"
function _dirFiles(path::String=".",rgxPat::Regex=Regex("Part_\\d{4}.bi4"))
    files = readdir(path)
    #Operation on dirFiles instantly
    filter!(x -> occursin(rgxPat, x), files)
    return files
end


readBi4Array(typ::Bi4DataType, Bi4Files::String) = readBi4Array(typ, [Bi4Files])
function readBi4Array(typ::Bi4DataType, Bi4Files::Vector{String} = _dirFiles())

    key    = Bi4DataTypeDict[typ].key
    offset = Bi4DataTypeDict[typ].offset
    T      = Bi4DataTypeDict[typ].type
    ncol   = Bi4DataTypeDict[typ].ncol
   
    j  = similar(Vector{Vector{T}},axes(Bi4Files))

    Base.Threads.@threads for i in eachindex(j)
        #println("$typ | Iteration: " * lpad(string(i), 4, "0") * "|Reading: " * Bi4Files[i])
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
    # Multiply with 4 (sizeof(Float32)) here since UInt8 size, times number of particles, times
    # times number of columns gives the correct indices in the rf array for
    # Float32, Int32 etc.
    did_s = nid_e + 1 + sizeof(Int32)
    did_e = did_s - 1 + sizeof(T) * n * ncol
    

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

#
function _searchValue(str2Search::Vector{UInt8},strNeedle::String,seekCounter::Int,OutputType::DataType,nOutputs::Int=1)
    # Add Control Characters
    possibleControlCharacters = ["\f";"\v";"\b";"\x02";"\x16";"\x17"] #Use Char('\f') to see UInt8 value! https://www.rapidtables.com/code/text/ascii-table.html

    strNeedle_CC = "";
    hitIndex     = -1;
    for pCC in possibleControlCharacters
        strNeedle_CC = "\0"*strNeedle*pCC
        hitIndex     = Base._searchindex(str2Search,codeunits(strNeedle_CC),seekCounter)
        if hitIndex != 0
            break
        end
    end
    
    if hitIndex == 0
        @printf "StringNeedle: %s, was not found in file.\n" strNeedle
        return nothing
    end

    loc_a = hitIndex+ncodeunits(strNeedle_CC)+3 #+3 instead of +4 since we added one char in front of strNeedle!
    loc_b = loc_a+(nOutputs*sizeof(OutputType)-1)
    range_ab = loc_a:loc_b

    valR      = reinterpret(OutputType,str2Search[range_ab])

    if nOutputs == 1
        return valR[1]
    else
        return valR
    end
end

# Cannot read text yet, not advised to use yet
function readBi4_Head_Config()
    Bi4Head = _dirFiles(Regex("Part_Head"))

    file    = Bi4Head[1]

    # Import a full bi4 file as Array{UInt8,1}
    ft = open(file, read = true)
    rf = read(ft)
    close(ft)

    searchVar =
        Dict("ViscoType"=>(1,UInt32),
             "ViscoValue"=>(1,Float32),
             "ViscoBoundFactor"=>(1,Float32),
             "Splitting"=>(1,UInt8),
             "Dp"=>(1,Float64),
             "H"=>(1,Float64),
             "B"=>(1,Float64),
             "RhopZero"=>(1,Float64),
             "MassBound"=>(1,Float64),
             "MassFluid"=>(1,Float64),
             "Gamma"=>(1,Float64),
             "Gravity"=>(3,Float32),
             "CasePosMin"=>(3,Float64),
             "CasePosMax"=>(3,Float64),
             "PeriXinc"=>(3,Float64),
             "PeriYinc"=>(3,Float64),
             "PeriZinc"=>(3,Float64),
             "Data2d"=>(1,UInt8),
             "Data2dPosY"=>(1,Float64),
             "Npiece"=>(1,UInt32),
             "FirstPart"=>(1,UInt32),
        )

    for (sV,sT) in searchVar
        n = first(sT)
        t = last(sT)
        if n == 1
            a =  _searchValue(rf,sV,1,t)
            @printf "%s: %2.6f |%s|\n" sV a t
        elseif n == 3
            a =  _searchValue(rf,sV,1,t,n)
            @printf "%s: %2.6f %2.6f %2.6f |%s|\n" sV a[1] a[2] a[3] t
        end
    end

end


function readBi4_Head(path::String=".")
    Bi4Head = _dirFiles(path,Regex("Part_Head"))

    file    = Bi4Head[1]

    # Import a full bi4 file as Array{UInt8,1}
    ft = open(file, read = true)
    rf = read(ft)
    close(ft)

    Needle        = "ITEM"
    SearchNeedle  = codeunits(Needle)
    
    Item_Locs     = Vector{Int64}()
    ind        = 0
    
    while true
        ind  = Base._searchindex(rf, SearchNeedle, ind+1)
        if ind == 0
            break
        end
        push!(Item_Locs,ind)
    end

    popfirst!(Item_Locs)          #Remove 1 info Item
    popfirst!(Item_Locs)          #Remove 2 info Item
    push!(Item_Locs,length(rf)+1) #Add last range

    Item_Ranges = Vector{UnitRange{Int64}}()
    for i = 1:length(Item_Locs)-1
        push!(Item_Ranges,Item_Locs[i]:(Item_Locs[i+1]-1))
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
        Count =  _searchValue(rf_,"Count",1,Int32)

        MkType = _searchValue(rf_,"MkType",1,Int32)

        # Have to skip "MkBlocks" syntax.. Actually fixed now, since using control characters in _searchValue!
        Mk = _searchValue(rf_,"Mk",1,Int32)

        ActualType_ = searchType(rf_)


        dct[ival] = OrderedDict("Type"=>ActualType_,"MkType"=>MkType,"Mk"=>Mk,"Begin"=>Bi4_IdCount, "Count"=>Count)

        Bi4_IdCount += Count;
    end

    return dct

end

# Not advised to use yet
function readBi4_Info()
    Bi4Info = _dirFiles(Regex("PartInfo"))

    file    = Bi4Info[1]

    # Import a full bi4 file as Array{UInt8,1}
    ft = open(file, read = true)
    rf = read(ft)
    close(ft)

    Needle        = "ITEM"
    SearchNeedle  = codeunits(Needle)
    
    Item_Locs     = Vector{Int64}()
    ind           = 0
    
    while true
        ind  = Base._searchindex(rf, SearchNeedle, ind+1)
        if ind == 0
            break
        end
        push!(Item_Locs,ind)
    end

    popfirst!(Item_Locs)          #Remove 1 info Item
    push!(Item_Locs,length(rf)+1) #Add last range

    Item_Ranges = Vector{UnitRange{Int64}}()
    for i = 1:length(Item_Locs)-1
        push!(Item_Ranges,Item_Locs[i]:(Item_Locs[i+1]-1))
    end

    dct = Vector{OrderedDict}(undef,length(Item_Ranges))
    for (ival,valRange) in enumerate(Item_Ranges)
        rf_  = rf[valRange]
        Npok =  _searchValue(rf_,"Npok",1,Int32)
        Nout =  _searchValue(rf_,"Nout",1,Int32)
        Nptotal =  _searchValue(rf_,"Nptotal",1,Int32)
        RunTime =  _searchValue(rf_,"RunTime",1,Float64)
        TimeSim =  _searchValue(rf_,"timesim",1,Float64)
        TimeStep = _searchValue(rf_,"TimeStep",1,Float64)
        Step = _searchValue(rf_,"Step",1,Int32) #Wrong results, since it finds "TimeStep too"..

        dct[ival] = OrderedDict("Npok"=>Npok,"Nout"=>Nout,"Nptotal"=>Nptotal,"RunTime"=>RunTime,"timesim"=>TimeSim,"TimeStep"=>TimeStep,"Step"=>Step)
    end

    return dct
end

end #PostSPH

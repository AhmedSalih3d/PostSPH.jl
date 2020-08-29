__precompile__()

module PostSPH

using Printf #To construct string message easily
using StaticArrays

#Reading of XML files
include("ReadXML.jl")

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
        ## Make two dicts,   one for the enum -> string and one for the data type
        const searchString = Dict(Points => "POINTS ",Idp=>"POINT_DATA ",Vel=>"Vel ",Rhop=>"Rhop ",Mass=>"Mass ",Press => "Press ",Vol => "Vol ",Ace => "Ace ",Vor => "Vor ",Typ=>"Type ",Mk => "Mk ")
        const catType = Dict(Points => Float32,Idp=>Int32,Vel=>Float32,Rhop=>Float32,Mass=>Float32,Press =>Float32,Vol=>Float32,Ace=>Float32,Vor=>Float32,Typ=>Int8,Mk=>Int8)

    #Purpose of this is to ensure faster read speed by finding approximate location
    #of wanted array (typ) in vtk file. The approximate location is found in "PosTyp" variable.
    function _readVtkPos(filename::String,typ::Cat)
        fd::IOStream = open(filename, read=true)
        readuntil(fd, searchString[typ])

        breakPos = 0
        if eof(fd) == true
            breakPos = 1
            PosTyp   = NaN
        else
            if typ != Points
                PosTyp   = position(fd)-div(position(fd),20)
            else
                PosTyp = 0
            end
        end
        return PosTyp,breakPos
    end

    ##Read VTK and chose option
    function _readVtk(filename::String, typ::Cat,PosTyp::Number)
        #Open the specific file
        fd::IOStream = open(filename, read=true)
        #Read until the specific enum string has been found. Ie. Mk => "Mk " etc.
        seek(fd,PosTyp)
        readuntil(fd, searchString[typ])
        #State number of colums depending on binary vtk file
        nCol::Int = if ( typ == Points ) 3 elseif (typ == Idp) 1 else parse(Int,readuntil(fd,' ')) end
        #Do something if special case?
        speCase::Bool = typ == Idp
        char = (speCase) ? '\n' : ' '
        nRow::Int = parse(Int,readuntil(fd,char))

        if ( speCase)
            readuntil(fd,'\n')
        end

        readuntil(fd,'\n')

        if (nCol > 1)
            dim = 2
            size = (nRow,nCol)
        else
            dim = 1
            size = nRow
        end
        #Preallocate an array depending on datatype and of chosen size
        #arrayVal::Array{catType[typ],dim} = zeros(catType[typ],size)
        arrayVal = Array{catType[typ],dim}(undef, size)

        _transferData(fd, arrayVal)

        #Close the open file
        close(fd)
        return arrayVal
    end

    #_transferData function used to read file depending on data type.

    function _transferData(fd::IOStream, arrayVal::Array{Float32,1})
        sz = size(arrayVal)
        for i = 1:sz[1]::Number
            @inbounds arrayVal[i] = ntoh(read(fd, Float32))
        end
    end

    function _transferData(fd::IOStream, arrayVal::Array{Float32,2})
        sz = size(arrayVal)
        @inbounds for i = 1:sz[1]::Number
            @inbounds for k = 1:sz[2]::Number
                @inbounds arrayVal[i,k] = ntoh(read(fd, Float32))
                      end
                  end
    end

    function _transferData(fd::IOStream, arrayVal::Array{Int32,1})
        sz = size(arrayVal)
        for i = 1:sz[1]::Number
            @inbounds arrayVal[i] = ntoh(read(fd, Int32))
        end
    end

    function _transferData(fd::IOStream, arrayVal::Array{Int8,1})
        sz = size(arrayVal)
        for i = 1:sz[1]::Number
            @inbounds arrayVal[i] = read(fd, Int8)
        end
    end

    ## Function to read multiplefiles into array
    #filename is a part of the wanted files - CASE SENSITIVE
    function readVtkArray(filename::String, typ::Cat)
        #Generate file names with dir command
        dirFiles  = readdir()
        filenames = filter!(s->occursin(filename, s),dirFiles) #Can't use r?
        nFilenames = size(filenames)[1]
        PosTyp,breakPos = _readVtkPos(filenames[1],typ)
        if breakPos == 1
            strMsg = @sprintf "%s was not found in .vtk file" typ
            println(strMsg)
            k = nothing
        else
            k = Vector{Array{catType[typ]}}(undef, nFilenames)
            Threads.@threads for i = 1:nFilenames::Number
                try
                    @inbounds k[i] = _readVtk(filenames[i], typ,PosTyp)
                catch
                    #Since DualSPHysics starts from 0000 - Test
                    println("Error in file number ",i-1)
                    @inbounds k[i] = NaN
                end
            end
        end
        return k
    end

###########################Extra Functionality#################################
    #Purpose is to only read number of particles in a simulation step. This
    #automatically uses "typ" = Points, since Points will always exist. It also
    #starts from PosTyp = 0, so this argument is not relevant anymore
    function _GetParticles(filename::String)
        typ = Points
        fd::IOStream = open(filename, read=true)
        #Read until the specific enum string has been found. Ie. Mk => "Mk " etc.
        readuntil(fd, searchString[typ])
        #State number of colums depending on binary vtk file
        nRow::Int = parse(Int,readuntil(fd,' '))
        return nRow
    end

    ## Read number of particles through each simulation step in Int64 Array
    function readVtkParticles(filename::String)
        #Generate file names with dir command
        dirFiles  = readdir()
        filenames = filter!(s->occursin(filename, s),dirFiles) #Can't use r?
        nFilenames = size(filenames)[1]
        k = Array{Int64}(undef, nFilenames)
Threads.@threads for i = 1:nFilenames::Number
                    @inbounds k[i] = _GetParticles(filenames[i])
                 end
        return k
    end

    #Function which reads all available variables in a SINGLE vtk file. "reset"
    #is set to "true" to ensure correct results, even though the order in "Cat",
    #should be true for all vtk files. Set to "false" if the order in "Cat" is right
    #and to get a 7.5 performance boost, ie. wanting to check 10000 files in a for loop.
    function readVtkVariables(filename::String,resetVar::Bool=true)
        k::Array{String} = []
        #Open the specific file
        fd::IOStream = open(filename, read=true)
        for i in instances(Cat)
            #Read until the specific enum string has been found. Ie. Mk => "Mk " etc.
            if resetVar == true
                seek(fd,0)
            end
            readuntil(fd, searchString[i])
            #Check if end of file (ie. variable not found)
            if eof(fd) != true
                k = push!(k,string(i))
            end
        end
        return k
    end

    #Function which returns available vtk files in current folder as a string array
    #Link: https://discourse.julialang.org/t/help-me-make-a-regex-instead-of-this/23878/7
    #Altered to fit coding convention in PostSPH. Able to insert own path/pattern if wanted.
    function readVtkNames(filenames::Array{String,1}=readdir(), pattern::Regex=r"([a-zA-Z]+).*\.vtk$")
        namelist = Vector{String}()  # this makes sure that the list is unique
        for filename in filenames
            m = match(pattern, filename)
            isnothing(m) && continue
            str = first(m.captures) #? captures for help, it just gets the "PartFluid", component
            !isempty(str) && push!(namelist, str)
        end
        return unique!(namelist)
    end

##Read from binary files directly
# Utilizes same "Cat"

## Make two dicts,   one for the enum -> string and one for the data type - Bi4 specific
const searchStringBi4 = Dict{Cat,String}(Idp => "ARRAY\x03\0\0\0Idp", Points => "ARRAY\x03\0\0\0Pos", Vel => "ARRAY\x03\0\0\0Vel", Rhop => "ARRAY\x04\0\0\0Rhop")
const catTypeBi4 = Dict{Cat,DataType}(Idp => Int32, Points => Float32, Vel => Float32, Rhop => Float32)
const catArrayBi4  = Dict{Cat,Int64}(Idp => 1, Points => 2, Vel => 2, Rhop => 1)
const catColBi4  = Dict{Cat,Int64}(Idp => 1, Points => 3, Vel => 3, Rhop => 1)

##Lists files in directory and only returns applicable files, ie. "Part_XXXX.bi4"
function _dirFiles()
    files = readdir()
    #Operation on dirFiles instantly
    filter!(x->occursin(r"Part_\d{4}.bi4",x),files)
    return files
end

function isc(char)
           char == "ARRAY\x03\0\0\0Idp"
end
##Function to determine correct location of array in file for first file and used
# as an approximate for all future files. Can be turned off.
# ALLOCATIONS IN _dirFiles() due to "readdir()" which is Julia base
function _Bi4Pos(typ::Cat,first_file::String)
    ft = open(first_file,read=true)
    @time readuntil(ft,searchStringBi4[typ])
    typPos = position(ft)-div(position(ft),20)
    close(ft)
    return typPos
end


##transferData is smart since we change the value of an array in place and
# skip a lot of unnecessary allocation steps. It has been tuned to accept
# matrices and vectors, where the type is automatically deduced using "eltype".
@inline _readel(ft::IOStream, typ) = read(ft, typ)
#@inline _readel(ft::IOStream, typ::Type{<:SVector{3,T}}) where T = SVector(_readel(ft, T), _readel(ft, T), _readel(ft, T))

function _transferDataBi4(ft::IOStream, arrayVal::AbstractVector)
    typ = eltype(arrayVal)

    @inbounds for i in eachindex(arrayVal)
        arrayVal[i] = _readel(ft, typ)
    end
end


##If nCol is not 1, then type is static array
_typeMaker(typ, nCol) = nCol == 1 ? typ : SVector{nCol, typ}
#If condition is true, set starting position to 0 of file.
_StartFromTop(condition,startPos) = condition == false ? startPos : 0

#Command for plotting
#scatter(getindex.(a[1], 1), getindex.(a[1], 3))
#readBi4Array(typ::Cat,Bi4Files::Array{String,1}=_dirFiles()) = readBi4Array(typ, false, Bi4Files)
readBi4Array(typ::Cat,Bi4Files::Array{String,1}) = readBi4Array(typ, false, Bi4Files)
readBi4Array(typ::Cat,SeekNull::Bool,Bi4Files::String) = readBi4Array(typ, false, [Bi4Files])
readBi4Array(typ::Cat,Bi4Files::String) = readBi4Array(typ, false, [Bi4Files])
function readBi4Array(typ::Cat,SeekNull::Bool=false,Bi4Files::Array{String,1}=_dirFiles())

    startPos = _Bi4Pos(typ,Bi4Files[1])
    #if true start from top, else do not, use startPos from _Bi4Pos
    _StartFromTop(SeekNull,startPos) = condition == false ? startPos : 0
    nBi4     = size(Bi4Files)[1]

    T  = _typeMaker(catTypeBi4[typ], catColBi4[typ])
    j  = fill(Array{T,1}(), nBi4) #Less allocs than Vector{Array{T}}(undef,nBi4)
    #j  = Vector{Array{T}}(undef,nBi4)

    Threads.@threads for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        seek(ft,startPos)
        readuntil(ft,searchStringBi4[typ])

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



end #PostSPH

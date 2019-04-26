__precompile__()

module PostSPH
export
    readVtkArray,
    Cat,
    ForceVtk
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


    function readVtkPos(filename::String,typ::Cat)
        fd::IOStream = open(filename, read=true)
        readuntil(fd, searchString[typ])
        if typ != Points
            PosTyp = position(fd)-div(position(fd),20)
        else
            PosTyp = 0
        end
        return PosTyp
    end
    ##Read VTK and chose option
    function readVtk(filename::String, typ::Cat,PosTyp::Number)
        #Open the specific file
        fd::IOStream = open(filename, read=true)
        #Read until the specific enum string has been found. Ie. Mk => "Mk " etc.
        seek(fd,PosTyp)
        readuntil(fd, searchString[typ])
        #State number of colums depending on binary vtk file
        nCol::UInt = if ( typ == Points ) 3 elseif (typ == Idp) 1 else parse(UInt,readuntil(fd,' ')) end
        #Do something if special case?
        speCase::Bool = typ == Idp
        char = (speCase) ? '\n' : ' '
        nRow::UInt = parse(UInt,readuntil(fd,char))

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
        arrayVal::Array{catType[typ],dim} = Array{catType[typ],dim}(undef, size)

        transferData(fd, arrayVal)

        #Close the open file
        close(fd)
        return arrayVal
    end

    #transferData function used to read file depending on type.

    function transferData(fd::IOStream, arrayVal::Array{Float32,1})
        sz = size(arrayVal)
        for i = 1:sz[1]::Number
            @inbounds arrayVal[i] = ntoh(read(fd, Float32))
        end
    end

    function transferData(fd::IOStream, arrayVal::Array{Float32,2})
        sz = size(arrayVal)
        @inbounds for i = 1:sz[1]::Number
            @inbounds for k = 1:sz[2]::Number
                @inbounds arrayVal[i,k] = ntoh(read(fd, Float32))
                      end
                  end
    end

    function transferData(fd::IOStream, arrayVal::Array{Int32,1})
        sz = size(arrayVal)
        for i = 1:sz[1]::Number
            @inbounds arrayVal[i] = ntoh(read(fd, Int32))
        end
    end

    function transferData(fd::IOStream, arrayVal::Array{Int8,1})
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
        k = Vector{Array{catType[typ]}}(undef, nFilenames)
        PosTyp = readVtkPos(filenames[1],typ)
        Threads.@threads for i = 1:nFilenames::Number
            try
                @inbounds k[i] = readVtk(filenames[i], typ,PosTyp)
            catch
                #Since DualSPHysics starts from 0000 - Test
                println("Error in file number ",i-1)
            end
        end

        return k
    end

    #These functions calculates force for either predefined arrays or filename using
    #mass of particle times acceleration of a particle. Use SplitVtk to get x y z
    #components and MagVtk to get magnitude of force later. Arrays have to be same
    #size which should always be true from data files, unless mistake in datafiles
    #x-force is given by Force = ForceVtk(filename) -->  Force[1][:,1], while
    #ForceMag is given by Force[2]
    function ForceVtk(filename::String)
        mass = readVtkArray(filename,Mass)
        ace  = readVtkArray(filename,Ace)
        n = size(mass)[1]
        ForceArray  = zeros(Float32,n,3)
        for i = 1:n
            @inbounds ForceArray[i,:] = sum(mass[i].*ace[i],dims=1)
        end
        ForceMag  = zeros(Float32,n)
        for i = 1:n
            @inbounds ForceMag[i] =  sqrt(ForceArray[i,1].^2+ForceArray[i,2].^2+ForceArray[i,3].^2)
        end
        return ForceArray,ForceMag
    end
end #PostSPH

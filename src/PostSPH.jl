__precompile__()

module PostSPH

using Printf #To construct string message easily

export
    readVtkArray,
    Cat,
    readVtkParticles,
    ForceVtk,
    MassVtk,
    readVtkNames,
    readVtkVariables
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
        arrayVal::Array{catType[typ],dim} = Array{catType[typ],dim}(undef, size)

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

    #Function which returns an array of mass of each single time step
    function MassVtk(filename::String)
        mass = readVtkArray(filename,Mass)
            if mass == nothing
                MassArray = nothing
            else
                n = size(mass)[1]
                MassArray  = Array{Float32,1}(undef,n)
                for i = 1:n
                    @inbounds MassArray[i] = sum(mass[i])
                end
            end
        return MassArray
    end

    #Arrays have to be same size which should always be true from data files,
    #unless mistake in datafiles x-force is given by
    #Force = ForceVtk(filename) -->  Force[1][:,1], while
    #ForceMag is given by Force[2]
    #Currently a bug exists which makes it so, if there is an error in one file
    #the whole terminal might crash..
    function ForceVtk(filename::String)
        mass = readVtkArray(filename,Mass)
        ace  = readVtkArray(filename,Ace)
            if mass == nothing || ace == nothing
                ForceArray = nothing
                ForceMag   = nothing
            else
                n = size(mass)[1]
                ForceArray  = Array{Float32,2}(undef,n,3)
                for i = 1:n
                    @inbounds ForceArray[i,:] = sum(mass[i].*ace[i],dims=1)
                end
                ForceMag  = Array{Float32,1}(undef,n)
                for i = 1:n
                    @inbounds ForceMag[i] =  sqrt(ForceArray[i,1].^2+ForceArray[i,2].^2+ForceArray[i,3].^2)
                end
            end
        return ForceArray,ForceMag
    end

    #Function to extract translation
    function FloatingVtkTranslation(filename::String)
        Pos = readVtkArray(filename,Points)
        n = length(Pos)
        k = size(Pos[1])[1] #Get number of particles in array
        x = Array{Float32,1}(undef,n)
        y = similar(x)
        z = similar(x)
         for i = 1:length(Pos)
             xyz = sum(Pos[i],dims=1)/k
             x[i] = xyz[1]
             y[i] = xyz[2]
             z[i] = xyz[3]
         end
         return x,y,z
     end

end #PostSPH

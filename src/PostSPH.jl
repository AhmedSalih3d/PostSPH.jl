__precompile__()

module PostSPH

using Printf #To construct string message easily

export
    readVtkArray,
    Cat,
    readVtkParticles,
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

    #_transferData function used to read file depending on type.

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

    #These functions calculates force for either predefined arrays or filename using
    #mass of particle times acceleration of a particle. Use SplitVtk to get x y z
    #components and MagVtk to get magnitude of force later. Arrays have to be same
    #size which should always be true from data files, unless mistake in datafiles
    #x-force is given by Force = ForceVtk(filename) -->  Force[1][:,1], while
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
end #PostSPH

###############################KERNELS##########################################
#h  = convert(Float32,0.002449) #m, smoothing length for PartFinal, River1
#x  = convert(Float32,0.2)
#\alphaD from Runout file
#αD = 92833.429688 #7/(4*pi*(h/2)^2) #Wendland Quintic 2D
#Abs since a particle on both sides of x should be counted in bool later
#q(ra::Float32,rb::Array{Float32,1},h::Float32) = abs.(ra.-rb)/h
##Wendland kernel as in Periodicity
#W(q) =  αD*(1 .- q/2).^4 .* (2 .* q .+ 1)
##Reading points data for particles and velocity
#pos = readVtkArray("PartFluid_",PostSPH.Points)
#vel = readVtkArray("PartFluid_",PostSPH.Vel)
##Initializing x-velocity array, Vab
#Vab = []
#for i = 1:length(pos)
#    qab = q(x,pos[i][:,1],h)
#    #Bool is to ensure on values between 0 and 2 are taken into account
#    bool = convert(Array{Int},0 .< qab .< 2)
#    qabF = q(x,pos[i][:,1],h).*bool
#    Wab = W(qabF)
#    push!(Vab,sum(Wab.*vel[i][:,1])/sum(Wab))
#end

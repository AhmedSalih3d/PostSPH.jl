#if typ == Idp || typ == Rhop
#    j[i] = zeros(catTypeBi4[typ], n)
#    _transferDataBi4(ft,j[i])
#elseif typ == Points || typ == Vel
#    j[i] = zeros(catTypeBi4[typ], (n,catColBi4[typ]))
#    _transferDataBi4(ft,j[i])
#end
##Constructs the array dimensions for the subarrays in the final vector, "j",
#see "readBi4Array"
function dimMaker(nRow::Int32,nCol::Int64)
    if nCol==1
        dim = (nRow,)
    elseif nCol==3
        dim = (nRow,nCol)
    end
    return dim
end

function _transferDataBi4(ft::IOStream, arrayVal::AbstractMatrix)
    typ = eltype(arrayVal)
    sz = size(arrayVal)
    if !eof(ft)
        for i = 1:sz[1]
            for k = 1:sz[2]
                @inbounds arrayVal[i,k] = read(ft, typ)
            end
        end
    end
end

function _transferDataBi4(ft::IOStream, arrayVal::AbstractVector)
    typ = eltype(arrayVal)
    sz = length(arrayVal)
    if !eof(ft)
        @inbounds for i in eachindex(arrayVal)
            arrayVal[i] = read(ft, typ)
        end
    end
end

function readBi4Array(typ::Cat,StartFromTop::Bool=false)
    Bi4Files = _dirFiles()

    if StartFromTop == false
        breakPos = _Bi4Pos(typ)
    else
        breakPos = 0
    end

    nBi4     = size(Bi4Files)[1]

    #j = Vector{Array{catTypeBi4[typ],catArrayBi4[typ]}}(undef, nBi4)
    j  = Vector{Array{}}(undef,nBi4)

    for i = 1:nBi4
        ft = open(Bi4Files[i],read=true)
        seek(ft,breakPos)
        readuntil(ft,searchStringBi4[typ])

            read(ft,Int64)
        n = read(ft,Int32)
            read(ft,Int32)
            println(n)
            dim = dimMaker(n,catColBi4[typ])
            j[i] = zeros(catTypeBi4[typ], dim)
            _transferDataBi4(ft,j[i])
        close(ft)
    end
    return j
end

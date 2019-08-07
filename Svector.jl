using StaticArrays

function mv(arr::Vector{T}, ::Val{N}) where {N,T}
       m=MVector{N,T}(undef)
       for i=1:N
       @inbounds m[i]=arr[i]
       end
       return SVector(m)
end

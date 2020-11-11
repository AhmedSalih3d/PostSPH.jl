using Distances
import Distances.colwise
using BenchmarkTools
usin

struct P{T}
    x::T
    y::T
    z::T
end

p_arr = fill(P{Float64}(rand(),rand(),rand()),100)

# Calculate distance between two vectors
x1 = [1.5,0,3]
y1 = [2.0,0,3]
@btime r1 = evaluate(Euclidean(), $x1, $y1)
# Result is "0.5" as exptected
#faster than sqrt(sum((x.-y).^2))

# Calculate colwise distance between two vectors
# Small test to understand if it works
# Has to be converted into "row" matrix
# Still "1D"
X1 = reshape(x1,(length(x1)),1)
Y1 = reshape(y1,(length(y1)),1)
@btime r2 = colwise(Euclidean(), X1, Y1)
# Testing with larger matrix - insert X1 and Y1 for easy check of accuracy
X1L = rand(3,100)
X1L[:,1] = X1
Y1L = rand(3,100)
Y1L[:,1] = Y1
@btime r2l = colwise(Euclidean(), $X1L, $Y1L)
# Test if it is faster to do manually
# This was not faster much more allocations. USING VIEWS FIXED IT SLIGHTLY
function calc_euclidean(X,Y)
    N = size(X,2)
    arr = Vector{Float64}(undef, N)
    @inbounds for i in eachindex(arr)
        arr[i] = @views evaluate(Euclidean(), view(X,:,i), view(Y,:,i))
    end

    return arr
end
@btime r2l_calc_euclidean = calc_euclidean($X1L,$Y1L)

N = 10^6
x2 = rand(N)
y2 = rand(N)
X2 = reshape(x2,(1,length(x2)))
Y2 = reshape(x2,(1,length(x2)))
@btime colwise(Euclidean(), $X2, $Y2)
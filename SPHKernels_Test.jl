using SPHKernels

# Kernel from here, I took Quintic (1 liner, love it like that)
#https://github.com/DualSPHysics/DualSPHysics/wiki/3.-SPH-formulation#31-smoothing-kernel

const αD1 = NaN
const αD2 = Float64(7/(4π))
const αD3 = Float64(21/(16π))

# Why is n_neighbours important?
# How would we know before hand, this can only be a rough estimate?
struct DSPH_Quintic{T} <: SPHKernel
    n_neighbours::Int64
    norm_1D::T
    norm_2D::T
    norm_3D::T
end

DSPH_Quintic(T::DataType=Float64, n_neighbours::Integer=216) = DSPH_Quintic{T}(n_neighbours, αD1, αD2, αD3)

@inline function SPHKernels.kernel_value_2D(kernel::DSPH_Quintic{T}, u::Real, h_inv::Real) where T

    # Option to disable if statement to speed up code?
    # Sometimes one can use "strict" neighbour finding in post-processing
    # which makes this condition obsolete.
    # After Thought: This way of writing is very fast anyways.
    if u < 2.0
        n = kernel.norm_2D * h_inv^2

        u_m1 = (1.0-u/2.0)
        u_m2 = (2.0*u+1.0)

        return (u_m1 * u_m1 * u_m1 * u_m1) * u_m2 * n |> T
    else
        return 0.0 |> T
    end

end

# I think that you potentially might have had a mathematical misunderstanding here
# as I have had my self in the past. The gradient of the kernel is not equal to the derivative 
# of the kernel. Because of the chain rule

# The derivative: dWdq (I would call this kernel_deriv)

# The gradient(W) = dWdq * (x - y)/|x-y| * (1/h) (I would call this kernel_grad)

# If you think I am wrong, please do challenge me on this point, I want to learn too! :-)

@inline function kernel_grad_2D(kernel::DSPH_Quintic{T}, u::Real, h_inv::Real,x::AbstractArray,y::AbstractArray,r::Real) where T
    

    n = kernel.norm_2D * h_inv^2


    # I am used to defining u as q
    # The derivative of the kernel, W(q), defined as dWdq is given in this case as:

    # W(q)    = (1-q/2)^4 * (2q + 1)

    # dWdq(q) = (5/8) * q * (q-2)^3


    if u < 2.0
        u_d1 = 5.0/8.0
        u_d2 = (u - 2.0)

        dWdq = u_d1 * u * (u_d2 * u_d2 * u_d2) |> T
    else
        return @. [0.0;0.0] |> T
    end

    # Calculate result for each direction

    x_diff = x[1] - y[1]
    y_diff = x[2] - y[2]

    # Notice that terms are common except "_diff"
    tmp = dWdq*r*h_inv

    Wgx = tmp*x_diff |> T
    Wgy = tmp*y_diff |> T

    # Return values

    return @. [Wgx;Wgy] |> T
end

# Play

k     = DSPH_Quintic()
u1    = 2.0
h     = 0.028284271247
h_inv = 1.0/h;

@time  v = kernel_value_2D(k, u1, h_inv)

x     = [2.0;2.00]
y     = [2.0;2.15]
r     = sum((x.-y).^2)
u2    =  r/h

Wgx,Wgy =  kernel_grad_2D(k, u2, h_inv,x,y,r)

# And to check consistency we flip x and y_diff
Wgx_f,Wgy_f =  kernel_grad_2D(k, u2, h_inv,y,x,r)

# Same result but opposite signs! This is more clearly shown when x[1] != y[1], disregard sign on "-0".

## BenchmarkTools

# @btime  v = kernel_value_2D($k, $u1, $h_inv)
# @btime  v = kernel_value_2D($k, $u1, $h_inv)
#   0.001 ns (0 allocations: 0 bytes)
# 0.0

# @btime Wgx,Wgy =  kernel_grad_2D($k, $u2, $h_inv,$x,$y,$r)
#   47.827 ns (2 allocations: 192 bytes)

# 2 allocs match allocating Wgx and Wgy, so seems like best possible case, but I don't understand how kernel_value can have absolute zero allocations
# Can you see an improvement in the grad function?
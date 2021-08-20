# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
rho   = rand(990:rand():1010,10^3)

function EOS1(rho;c0=100,gamma=7,rho0=1000)
    b=(c0^2*rho0)/gamma

    P     = b*((rho/rho0).^7 .- 1)
end

function EOS1SIMILAR(rho;c0=100,gamma=7,rho0=1000)
    b=(c0^2*rho0)/gamma

    P     = similar(rho)

    P     .= b*((rho/rho0).^7 .- 1)
end

function EOS2(rho;c0=100,gamma=7,rho0=1000)
    b=(c0^2*rho0)/gamma

    rho0_denominator = 1/(rho0^7);

    P     = b*(rho.^7 * rho0_denominator .- 1)
end

function EOS3(rho;c0=100,gamma=7,rho0=1000)
    b=(c0^2*rho0)/gamma

    rho0_denominator = 1/(rho0^7);

    P     = rho.^7 * rho0_denominator * b .- b
end

bEOS1  = @benchmark EOS1($rho)
println("Minimum Time EOS1:\n",bEOS1)
bEOS2  = @benchmark EOS2($rho)
println("Minimum Time EOS2:\n",bEOS2)
bEOS3  = @benchmark EOS3($rho)
println("Minimum Time EOS3:\n",bEOS3)

bEOS1AVX  = @benchmark @avx EOS1($rho)
println("Minimum Time EOS1AVX:\n",bEOS1AVX)
bEOS1SIMILAR  = @benchmark EOS1SIMILAR($rho)
println("Minimum Time EOS1SIMILAR:\n",bEOS1SIMILAR)
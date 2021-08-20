# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
rho   = rand(990:rand():1010,10^3)
P     = similar(rho)

function EOSSimple(rho;c0=100.0,rho0=1000.0)
    return  c0^2*(rho-rho0)
end


function EOS(rho;c0=100.0,gamma=7,rho0=1000.0)
    b=(c0^2*rho0)/gamma
    return  b*((rho/rho0)^gamma - 1.0)
end

function evaluateFunc!(func,array,arr)
    @tturbo @. arr = func(array)
end

@benchmark evaluateFunc!($EOSSimple,$rho,$P)
@benchmark evaluateFunc!($EOS,$rho,$P)

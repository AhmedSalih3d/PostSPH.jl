# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
using NearestNeighbors
using LinearAlgebra
#using Plots

N     = 10^3 
q     = rand(0:rand():2,N)
rho   = rand(990:rand():1010,N)

x     = similar(rho)
P     = similar(rho)

function EOSSimple(rho;c0=100.0,rho0=1000.0)
    return  c0^2*(rho-rho0)
end


function EOS(rho;c0=100.0,gamma=7,rho0=1000.0)
    b=(c0^2*rho0)/gamma
    return  b*((rho/rho0)^gamma - 1.0)
end

# CubicSpline not easy to implement.
function Wendland(q;aD=1.0)
    
    return aD*(1.0-q/2.0)^4 * (2.0*q+1.0)
end

function evaluateFunc!(func,array,arr)
    @tturbo @. arr = func(array)
end

@benchmark evaluateFunc!($EOSSimple,$rho,$P)
@benchmark evaluateFunc!($EOS,$rho,$P)

## EXAMPLE
using PostSPH

# Load all data in
cd(raw"D:\DualSPHysics_v5.0\examples\main\03_MovingSquare\CaseMovingSquare_out\data")
rhop_array = readBi4Array(PostSPH.Rhop)
pos_array  = readBi4Array(PostSPH.Points)
vel_array  = readBi4Array(PostSPH.Vel)
idp_array  = readBi4Array(PostSPH.Idp)
##

for it = 1:251
    xx = pos_array[it][1:3:end]
    yy = pos_array[it][2:3:end]
    zz = pos_array[it][3:3:end]

    data = hcat(xx,yy,zz)'
    H    = 0.028284271247
    TOH  = 2.0*H

    balltree = BallTree(data)
    idxs     = inrange(balltree, data, TOH, true)

    r_mag    = map(x->norm.(eachcol(data[:,x] .- data[:,1])),idxs)

    for i    = 1:length(idxs)
        r_mag[i] = norm.(eachcol( data[:,i] .- data[:,idxs[i]]))
    end

    q        = r_mag./H;
    clamp!.(q,0.0,2.0)

    aD       = (7.0)/(4.0*π*H^2) #2d

    Wab      = map(x -> sum(Wendland.(x,aD=aD)),q)

    SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array[it],
    Idp    = Wab,
    Vel    = vel_array[it],
    Rhop   = rhop_array[it])

    #Use 3d glyph filter in Paraview with 2d glyphs!
    PostSPH.SaveVTK.write_vtp("SimDataYay_"*lpad(string(it),4,"0"),SimData101)
end
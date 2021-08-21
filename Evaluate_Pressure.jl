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

function WendlandDerivative(q;aD=1.0)
    
    return aD*((5.0/8.0)*q*(q-2.0)^3)
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


H    = 0.028284271247
TOH  = 2.0*H
aD   = (7.0)/(4.0*π*H^2) #2d

ϵ    = 1e-6;
mb   = 0.4;

function WandWg(it,pos_array,rhop_array,H)
    xx = pos_array[it][1:3:end]
    yy = pos_array[it][2:3:end]
    zz = pos_array[it][3:3:end]

    data = hcat(xx,yy,zz)'
    
    np = length(xx)

    balltree = BallTree(data)
    idxs     = inrange(balltree, data, TOH, true)

    r_mag    = Array{Array{Float32,1},1}(undef,np)
    
    @time  for i    = 1:np
        r_mag[i] = norm.(eachcol( data[:,i] .- data[:,idxs[i]]))
    end

    q        = r_mag./H;
    clamp!.(q,0.0,2.0)

    Wab      = map(x -> sum(Wendland.(x,aD=aD)),q)

    WabM     = map(x -> sum(Wendland.(x,aD=aD)),q)

    WabM = similar(Wab)
    for i = 1:np
        WabM[i] =  sum(Wendland.(q[i],aD=aD)./sum(mb./rhop_array[it][idxs[i]] .* Wendland.(q[i],aD=aD)))
    end

    Wgx      = Array{Float32,1}(undef,np)
    Wgy      = zeros(Float32,np);#Array{Float32,1}(undef,np)
    Wgz      = Array{Float32,1}(undef,np)
    for i = 1:np
        Wgx[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[1,i] .- data[1,idxs[i]]) ./ (r_mag[i] * H .+ ϵ ))
        #Wgy[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[2,i] .- data[2,idxs[i]]) ./ (r_mag[i] .+ H ))
        Wgz[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[3,i] .- data[3,idxs[i]]) ./ (r_mag[i] * H .+ ϵ ))
    end
    
    Wg      = zeros(Float32,3*np)
    
    Wg[1:3:end] = Wgx
    Wg[2:3:end] = Wgy
    Wg[3:3:end] = Wgz   

    SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array[it],
    Idp    = Wab,
    Vel    = Wg,
    Rhop   = WabM*mb)

    #Use 3d glyph filter in Paraview with 2d glyphs!
    PostSPH.SaveVTK.write_vtp("SimDataYay_"*lpad(string(it),4,"0"),SimData101)
end

for it = 1:length(rhop_array)
    @time WandWg(it,pos_array,rhop_array,H)
end
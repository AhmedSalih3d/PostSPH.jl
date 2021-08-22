# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
using NearestNeighbors
using LinearAlgebra
using Base.Threads
using Distances
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
id_range   = (1501,4101)
rhop_array = readBi4Array(PostSPH.Rhop,id_range)
pos_array  = readBi4Array(PostSPH.Points,id_range)
vel_array  = readBi4Array(PostSPH.Vel,id_range)
idp_array  = readBi4Array(PostSPH.Idp,id_range)
##
Npok       = PostSPH.readBi4_CurrentTotalParticles()
np_max     = maximum(Npok)

H    = 0.028284271247
TOH  = 2.0*H
aD   = (7.0)/(4.0*π*H^2) #2d

ϵ    = 1e-6;
mb   = 0.4;
d1 = Dict{String,Array}()

function WandWg(it,pos_array,rhop_array,H)
    xx = @view pos_array[it][1:3:end]
    yy = @view pos_array[it][2:3:end]
    zz = @view pos_array[it][3:3:end]

    data = hcat(xx,yy,zz)'
    
    np = length(xx)

    balltree = BallTree(data)
    idxs     = inrange(balltree, data, TOH, true)

    r_mag    = Array{Array{Float32,1},1}(undef,np)
    
    @time @threads for i    = 1:np
        #r_mag[i] = norm.(eachcol( data[:,i] .- data[:,idxs[i]]))
        r_mag[i] = colwise(Euclidean(), data[:,i],@view data[:,idxs[i]])
    end

    q        = r_mag./H;
    clamp!.(q,0.0,2.0)

    Wab      = map(x -> sum(Wendland.(x,aD=aD)),q)

    WabM     = map(x -> sum(Wendland.(x,aD=aD)),q)

    WabM = similar(Wab)
    @time @threads for i = 1:np
        Wab_tmp = Wendland.(q[i],aD=aD)
        denom   =  mb*sum(1 ./ rhop_array[it][idxs[i]] .* Wab_tmp)
        WabM[i] =  sum(Wab_tmp./denom)
    end

    rhoM     = WabM*mb;

    Wgx      = Array{Float32,1}(undef,np)
    Wgy      = zeros(Float32,np);#Array{Float32,1}(undef,np)
    Wgz      = Array{Float32,1}(undef,np)
    @threads for i = 1:np
        Wgx[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[1,i] .- data[1,idxs[i]]) ./ (r_mag[i] * H .+ ϵ ))
        #Wgy[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[2,i] .- data[2,idxs[i]]) ./ (r_mag[i] .+ H ))
        Wgz[i] = sum(WendlandDerivative.(q[i],aD=aD) .* (data[3,i] .- data[3,idxs[i]]) ./ (r_mag[i] * H .+ ϵ ))
    end
    
    Wg      = zeros(Float32,3*np)
    
    Wg[1:3:end] = Wgx
    Wg[2:3:end] = Wgy
    Wg[3:3:end] = Wgz   

    SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array[it],
    Idp    = idp_array[it], #Wab,
    Vel    = vel_array[it], #Wg
    Rhop   = rhop_array[it])#rhoM)

    d1["Wab"]  = Wab
    d1["Wg"]   = Wg
    d1["rhoM"] = rhoM

    #Use 3d glyph filter in Paraview with 2d glyphs!
    PostSPH.SaveVTK.write_vtp("SimDataYay_"*lpad(string(it),4,"0"),SimData101,extra_arrays=d1)
end

for it = 1:length(rhop_array)
    @time WandWg(it,pos_array,rhop_array,H)
end
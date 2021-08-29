# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
using NearestNeighbors
using LinearAlgebra
using Base.Threads
using Distances
using ProfileView
#using Plots

# N     = 10^3 
# q     = rand(0:rand():2,N)
# rho   = rand(990:rand():1010,N)

# x     = similar(rho)
# P     = similar(rho)

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

# @benchmark evaluateFunc!($EOSSimple,$rho,$P)
# @benchmark evaluateFunc!($EOS,$rho,$P)

## EXAMPLE
using PostSPH

# Load all data in
cd(raw"D:\DualSPHysics_v5.0\examples\main\03_MovingSquare\CaseMovingSquare_out\data")
id_range   = ()#(1501,4101)
rhop_array = readBi4Array(PostSPH.Rhop)
pos_array  = readBi4Array(PostSPH.Points)
vel_array  = readBi4Array(PostSPH.Vel)
idp_array  = readBi4Array(PostSPH.Idp)
##
Npok       = PostSPH.readBi4_CurrentTotalParticles()
np_max     = maximum(Npok)

H    = 0.028284271247
TOH  = 2.0*H
aD   = (7.0)/(4.0*π*H^2) #2d
ϵ    = 1e-6;
mb   = 0.4;
d1 = Dict{String,Array}()

function constructData(pos_array,it)
    xx = @view pos_array[it][1:3:end]
    yy = @view pos_array[it][2:3:end]
    zz = @view pos_array[it][3:3:end]
    data = hcat(xx,yy,zz)'

    return data
end

function FindNeighbours(data)
    balltree = BallTree(data)
    idxs     = inrange(balltree, data, TOH, false) #can be set to false

    return idxs
end

function CalculateRmag(data,idxs)
    np = length(idxs)

    r_mag    = Array{Array{Float32,1},1}(undef,np)
    
    @threads for i    = 1:np
        #r_mag[i] = norm.(eachcol( data[:,i] .- data[:,idxs[i]]))
        r_mag[i] = colwise(Euclidean(), data[:,i],@view data[:,idxs[i]])
    end

    return r_mag
end

function CalculateQClamp(r_mag,H)
    q        = r_mag./H;
    clamp!.(q,0.0,2.0)

    return q
end

function CalculateShephardFilterMass(Wab,it,q,idxs)

    np    = length(idxs)

    WabM  = similar(Wab)
    @threads for i = 1:np
        Wab_tmp =  Wendland.(q[i],aD=aD)
        denom   =  mb*sum(1 ./ rhop_array[it][idxs[i]] .* Wab_tmp)
        WabM[i] =  sum(Wab_tmp./denom)
    end
    rhoM     = WabM*mb;

    return rhoM
end

function CalculateXIJ(data,idxs,r_mag)
    np = length(idxs)

    function tt(v, A)
        r = similar(A)
        @inbounds for j = 1:size(A,2) 
            @simd for i = 1:size(A,1) 
                r[i,j] = v[j] * A[i,j]
            end
        end 
        r
    end 

    xij      = Array{Array{Float32,1},1}(undef,np)
    yij      = similar(xij)
    zij      = similar(xij)
    @threads for i = 1:np
        denom            = 1 ./ (r_mag[i]*H.+ϵ)
        neighbour_ids    = idxs[i]
        p  = data[:,i]
        XIJ  = tt(denom,(p .- @view data[:,neighbour_ids]))

        xij[i] = XIJ[1,:]
        yij[i] = XIJ[2,:]
        zij[i] = XIJ[3,:]
    end

    return xij,yij,zij
end

function CalculateGradient(q,xij,zij)

    np       = length(q)

    Wg       = zeros(Float32,3*np)

    @threads for i = 1:np
        dWdq   = WendlandDerivative.(q[i],aD=aD)
        offset = i+(2*(i-1))
        Wg[offset]   = sum(dWdq .*  xij[i])
        #Wgy[i] = sum(WendlandDerivative.(q[i],aD=aD) .*  yij[i])
        Wg[offset+2] = sum(dWdq .*  zij[i])
    end
    
    return Wg
end

function WandWg(it,pos_array,rhop_array,H)

    data  = constructData(pos_array,it)

    idxs  = FindNeighbours(data)

    r_mag = CalculateRmag(data,idxs)

    q     = CalculateQClamp(r_mag,H)

    Wab   = map(x -> sum(Wendland.(x,aD=aD)),q)

    rhoM  = CalculateShephardFilterMass(Wab,it,q,idxs)

    xij,_,zij = CalculateXIJ(data,idxs,r_mag)

    
    Wg = CalculateGradient(q,xij,zij);

    SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array[it],
    Idp    = idp_array[it], 
    Vel    = vel_array[it], 
    Rhop   = rhop_array[it])

    d1["Wab"]  = Wab
    d1["Wg"]   = Wg
    d1["rhoM"] = rhoM

    #Use 3d glyph filter in Paraview with 2d glyphs!
    PostSPH.SaveVTK.write_vtp("SimDataYay_"*lpad(string(it),4,"0"),SimData101,extra_arrays=d1)
end

function WandWgOuter(pos_array,rhop_array,H)
    p   = pos_array;
    r   = rhop_array
    for it = 1:length(rhop_array)
        println(it)
        @time WandWg(it,p,r,H)
    end
end

#@profview

WandWgOuter(pos_array,rhop_array,H)


# Testing fastest implementation of evaluating scalar function and passing a scalar function to others
using BenchmarkTools
using LoopVectorization
using NearestNeighbors
using LinearAlgebra
using Base.Threads
using Distances
using ProfileView
using StaticArrays
using CellListMap
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
function Wendland(q,aD)
    
    return aD*(1.0-q/2.0)^4 * (2.0*q+1.0)
end

function WendlandDerivative(q,aD)
    
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

const H    = 0.028284271247
const TOH  = 2.0*H
const aD   = (7.0)/(4.0*π*H^2) #2d
const ϵ    = 1e-6;
const mb   = 0.4;
d1 = Dict{String,Array}()

function constructData(pos_array)
    xx = @view pos_array[1:3:end]
    yy = @view pos_array[2:3:end]
    zz = @view pos_array[3:3:end]
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

function CalculateShephardFilterMass(Wab,rhop_array,q,idxs)

    np    = length(idxs)

    WabM  = similar(Wab)
    @threads for i = 1:np
        Wab_tmp =  Wendland.(q[i],aD)
        denom   =  mb*sum(1 ./ rhop_array[idxs[i]] .* Wab_tmp)
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
        dWdq   = WendlandDerivative.(q[i],aD)
        offset = i+(2*(i-1))
        Wg[offset]   = sum(dWdq .*  xij[i])
        #Wgy[i] = sum(WendlandDerivative.(q[i],aD=aD) .*  yij[i])
        Wg[offset+2] = sum(dWdq .*  zij[i])
    end
    
    return Wg
end

function WandWg(it,pos_array,rhop_array,idp_array,vel_array,H)

    data  = constructData(pos_array)

    idxs  = FindNeighbours(data)

    r_mag = CalculateRmag(data,idxs)

    q     = CalculateQClamp(r_mag,H)

    Wab   = map(x -> sum(Wendland.(x,aD)),q)

    rhoM  = CalculateShephardFilterMass(Wab,rhop_array,q,idxs)

    xij,_,zij = CalculateXIJ(data,idxs,r_mag)

    
    Wg = CalculateGradient(q,xij,zij);

    SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array,
    Idp    = idp_array, 
    Vel    = vel_array, 
    Rhop   = rhop_array)

    d1["Wab"]  = Wab
    d1["Wg"]   = Wg
    d1["rhoM"] = rhoM

    #Use 3d glyph filter in Paraview with 2d glyphs!
    PostSPH.SaveVTK.write_vtp("SimDataYay_"*lpad(string(it),4,"0"),SimData101,extra_arrays=d1)
end

function WandWgOuter(pos_array,rhop_array,H)
    p   = pos_array;
    r   = rhop_array;
    id  = idp_array;
    v   = vel_array;
    for it = 1:length(rhop_array)
        println(it)
        @time WandWg(it,p[it],r[it],id[it],v[it],H)
    end
end

#@profview

WandWgOuter(pos_array,rhop_array,H)


###  
using CellListMap
using StaticArrays

# System properties
N = 100_000
sides = [250,250,250]
cutoff = 10

# Particle positions
x = [ sides .* rand(3) for i in 1:N ]

# Initialize linked lists and box structures
box = Box(limits(x),cutoff)
cl = CellList(x,box)

mass = rand(N)

# Function to be evalulated for each pair: build distance histogram
function calc_forces!(x,y,i,j,d2,mass,forces)
    G = 9.8*mass[i]*mass[j]/d2
    d = sqrt(d2)
    df = (G/d)*(x - y)
    forces[i] = forces[i] - df
    forces[j] = forces[j] + df
    return forces
end

# Initialize and preallocate forces
forces = [ zeros(SVector{3,Float64}) for i in 1:N ]

# Run pairwise computation
map_pairwise!(
    (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
    forces,box,cl
)

## Own play
pos = deepcopy(pos_array[1])
deleteat!(pos,2:3:length(pos))
x   = reinterpret(SVector{2, eltype(pos)}, pos)

N   = length(x)


box = Box(limits(x),TOH)
cl = CellList(x,box)

# Function to be evalulated for each pair: build distance histogram
function calc_q!(x,y,i,j,d2,H,aD,q_arr)
    d = sqrt(d2)
    q_arr[i] += Wendland(d/H,aD)
    q_arr[j] += Wendland(d/H,aD)
    return q_arr
end

function calc_grad!(x,y,i,j,d2,H,aD,idx,g_arr)
    d   = sqrt(d2)
    q   = d/H
    tmp = (WendlandDerivative(q,aD)/(d+1e-6 ))*(1/H)
    g_arr[i] += tmp*(x[idx]-y[idx])
    g_arr[j] += tmp*(y[idx]-x[idx])
    return g_arr
end

# Initialize and preallocate forces
q_arr = zeros(Float32,N) #.+ Wendland(0.0,aD)#[ zeros(SVector{1,Float64}) for i in 1:N ]


# Run pairwise computation
map_pairwise!(
    (x,y,i,j,d2,q_arr) -> calc_q!(x,y,i,j,d2,TOH*0.5,aD,q_arr),
    q_arr,box,cl
)

@time q_arr .= q_arr .+ Wendland(0.0,aD) #Since own self is not found



gx = zeros(Float32,N)
gz = zeros(Float32,N)
idx = 1;
idz = 2;
map_pairwise!(
    (x,y,i,j,d2,gx) -> calc_grad!(x,y,i,j,d2,TOH*0.5,aD,idx,gx),
    gx,box,cl
)

map_pairwise!(
    (x,y,i,j,d2,gz) -> calc_grad!(x,y,i,j,d2,TOH*0.5,aD,idz,gz),
    gz,box,cl
)
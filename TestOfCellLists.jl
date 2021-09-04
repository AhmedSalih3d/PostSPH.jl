using CellListMap
using StaticArrays

# Define Functions for Calculating
function Wendland(q,aD)
    
    return aD*(1.0-q/2.0)^4 * (2.0*q+1.0)
end

function WendlandDerivative(q,aD)
    
    return aD*((5.0/8.0)*q*(q-2.0)^3)
end

# System properties
#N = 100_000
#sides = [50,50,50]
#cutoff = 2

# Particle positions
#x = [ sides .* rand(3) for i in 1:N ]




const aD   = Float32(1/sqrt(π))


function calc_grad!(x,y,i,j,d2,H,aD,idx,g_arr)
    d   = sqrt(d2)
    q   = d/H
    tmp = (WendlandDerivative(q,aD)/(d+1e-6 ))*(1/H)
    g_arr[i] += tmp*(x[idx]-y[idx])
    g_arr[j] += tmp*(y[idx]-x[idx])
    return g_arr
end

# Initialize and preallocate forces
q_arr = zeros(Float32,maximum(Npok))

# Function to be evalulated for each pair: build distance histogram
function calc_Wab!(x,y,i,j,d2,H,aD,q_arr)
    d = sqrt(d2)
    q_arr[i] += Wendland(d/H,aD)
    q_arr[j] += Wendland(d/H,aD)
    return q_arr
end

# Play
pos  = pos_array[1]
deleteat!(pos,2:3:length(pos))
x = reinterpret(SVector{2, eltype(pos)}, pos)

# Initialize linked lists and box structures
box = Box(limits(x),cutoff)
cl = CellList(x,box)
aux = CellListMap.AuxThreaded(cl)
for iter = 1:2
    pos  = pos_array[iter]
    deleteat!(pos,2:3:length(pos))
    x = reinterpret(SVector{2, eltype(pos)}, pos)

    box = Box(limits(x),cutoff)
@time    cl = UpdateCellList!(x,box,cl,aux) 
map_pairwise!(
        (x,y,i,j,d2,q_arr) -> calc_Wab!(x,y,i,j,d2,H,aD,q_arr),
        q_arr,box,cl
    );
end

q_arr .= q_arr .+ Wendland(0.0,aD) #Since own self is not found

# Calculate Gradients
gx = zeros(Float32,N)
gz = zeros(Float32,N)
idx = 1;
idz = 2;
for _ = 1:10
map_pairwise!(
    (x,y,i,j,d2,gx) -> calc_grad!(x,y,i,j,d2,cutoff/2,aD,idx,gx),
    gx,box,cl
)
end

map_pairwise!(
    (x,y,i,j,d2,gz) -> calc_grad!(x,y,i,j,d2,cutoff/2,aD,idz,gz),
    gz,box,cl
)

# Initialize cell lists with initial coordinates
cl = CellList(x,box)
# Allocate auxiliary arrays for threaded cell list construction
aux = CellListMap.AuxThreaded(cl)
for i in 1:nsteps
    x = ... # new coordinates
    box = Box(sides,cutoff) # perhaps the box has changed
    cl = UpdateCellList!(x,box,cl,aux) 
end
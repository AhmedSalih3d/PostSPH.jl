
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
for i = 1:10
@time map_pairwise!(
    (x,y,i,j,d2,forces) -> calc_forces!(x,y,i,j,d2,mass,forces),
    forces,box,cl
);
end

using Base.Threads
forces = [ zeros(SVector{3,Float64}) for i in 1:N ]
forces_threaded = [ deepcopy(forces) for i in 1:nthreads() ]

for i = 1:10
@time map_pairwise!(
   (x,y,i,j,d2, forces, output_threaded=forces_threaded) -> calc_forces!(x,y,i,j,d2,mass,forces),
   forces,box,cl
);
end
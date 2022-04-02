module SaveVTK

using WriteVTK
using Parameters

export write_vtp

"""
    Struct to hold simulation data
"""
@with_kw mutable struct SimData{T,N}
    Points::Vector{T}
    Idp::Vector{N}
    Vel::Vector{T}
    Rhop::Vector{T}
end

"""
    write_vtp(filename,points,attribute)

Saves a polydata file (.vtp) given a filename, a 2xN point array and and 1xN attribute array.
"""
function write_vtp(filename::String, sim_arr::SimData,act_id)

    x = @view sim_arr.Points[1:3:end][act_id]
    y = @view sim_arr.Points[2:3:end][act_id]
    z = @view sim_arr.Points[3:3:end][act_id]

    points = hcat(x, y, z)'

    polys = empty(MeshCell{WriteVTK.PolyData.Polys,UnitRange{Int64}}[])
    verts = empty(MeshCell{WriteVTK.PolyData.Verts,UnitRange{Int64}}[])

    # Note: the order of verts, lines, polys and strips is not important.
    # One doesn't even need to pass all of them.
    all_cells = (verts, polys)

    #allocations here, not my library
    vtk = vtk_grid(filename, points, all_cells..., compress = true, append = false)


    if !isempty(sim_arr.Idp)
        vtk["Idp"] = sim_arr.Idp[act_id]
    end
    if !isempty(sim_arr.Vel)
        vx = sim_arr.Vel[1:3:end][act_id]
        vy = sim_arr.Vel[2:3:end][act_id]
        vz = sim_arr.Vel[3:3:end][act_id]
        velocity = hcat(vx, vy, vz)'
        vtk["Velocity"] = velocity
    end
    if !isempty(sim_arr.Rhop)
        vtk["Rhop"] = sim_arr.Rhop[act_id]
    end

    vtk_save(vtk)
end

end

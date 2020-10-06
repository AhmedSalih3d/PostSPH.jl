module SaveVTK

    using StaticArrays
    using WriteVTK
    using Random
    using Parameters

    export
        write_vtp

    @with_kw mutable struct SimData
        Points::Array{Float32,2} = Array{Float32}(undef, 0, 0)
        Idp::Array{Int32,1} = Array{Int32,1}()
        Vel::Array{Float32,2} = Array{Float32}(undef, 0, 0)
        Rhop::Array{Float32,1} = Array{Float32,1}()
    end

    """
        write_vtp(filename,points,attribute)

    Saves a polydata file (.vtp) given a filename, a 2xN point array and and 1xN attribute array.
    """
    function write_vtp(filename::String, sim_arr::SimData)

        points   = sim_arr.Points

        polys = empty(MeshCell{WriteVTK.PolyData.Polys,UnitRange{Int64}}[])
        verts = empty(MeshCell{WriteVTK.PolyData.Verts,UnitRange{Int64}}[])

        # Note: the order of verts, lines, polys and strips is not important.
        # One doesn't even need to pass all of them.
        all_cells = (verts, polys)

        #allocations here, not my library
        vtk  = vtk_grid(filename, points, all_cells..., compress=true, append=false)


        if !isempty(sim_arr.Idp)
            vtk["Idp"]      = sim_arr.Idp
        end
        if !isempty(sim_arr.Vel)
            velocity = sim_arr.Vel
            vtk["Velocity"] = velocity
        end
        if !isempty(sim_arr.Rhop)
            vtk["Rhop"]     = sim_arr.Rhop
        end

        vtk_save(vtk)
    end

end

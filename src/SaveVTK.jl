module SaveVTK

    using StaticArrays
    using WriteVTK
    using Random
    using Parameters

    export
        write_vtp

    """
        Struct to hold simulation data
    """
    @with_kw mutable struct SimData
        Points#::Array{Float32,1} = Array{Float32}(undef, 0)
        Idp#::Array{Int32,1} = Array{Int32,1}()
        Vel#::Array{Float32,1} = Array{Float32}(undef, 0)
        Rhop#::Array{Float32,1} = Array{Float32,1}()
    end

    """
        write_vtp(filename,points,attribute)

    Saves a polydata file (.vtp) given a filename, a 2xN point array and and 1xN attribute array.
    """
    function write_vtp(filename::String, sim_arr::SimData;extra_arrays=Dict())

        

        x        = sim_arr.Points[1:3:end]
        y        = sim_arr.Points[2:3:end]
        z        = sim_arr.Points[3:3:end]

        nx       = length(x)


        points   = hcat(x,y,z)'

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
            vx       = sim_arr.Vel[1:3:end]
            vy       = sim_arr.Vel[2:3:end]
            vz       = sim_arr.Vel[3:3:end]
            velocity = hcat(vx,vy,vz)'
            vtk["Velocity"] = velocity
        end
        if !isempty(sim_arr.Rhop)
            vtk["Rhop"]     = sim_arr.Rhop
        end
        if !isempty(extra_arrays)
            for (key,value) in extra_arrays
                if length(value) > nx
                    valx = value[1:3:end]
                    valy = value[2:3:end]
                    valz = value[3:3:end]
                    val  = hcat(valx,valy,valz)'
                    vtk[key] = val
                else
                    vtk[key] = value
                end
            end
        end

        vtk_save(vtk)
    end

end

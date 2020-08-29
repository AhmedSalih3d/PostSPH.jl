module SaveVTK

    using StaticArrays
    using WriteVTK
    using Random


    export
        write_vtp

    """
    	 @Name(arg)
    Return the name of a variable, used for naming data fields in vtk files
    the same as the input Julia data array
    """
    #Does not work as expected
    macro Name(arg)
       string(arg)
    end


    export write_vtp, Name

    """
        write_vtp(filename,points,attribute)

    Saves a polydata file (.vtp) given a filename, a 2xN point array and and 1xN attribute array.
    """
    function write_vtp(filename::String, points, attribute) where T<:Number

        nr = size(points)[1]
        nc = size(points[1])[1]

        points = collect(Iterators.flatten(Iterators.flatten(points)))

        points = reshape(points,(nc, nr))

        polys = empty(MeshCell{WriteVTK.PolyData.Polys,UnitRange{Int64}}[])
        verts = empty(MeshCell{WriteVTK.PolyData.Verts,UnitRange{Int64}}[])

        # Note: the order of verts, lines, polys and strips is not important.
        # One doesn't even need to pass all of them.
        all_cells = (verts, polys)

        fname = filename
        vtk  = vtk_grid(fname, points, all_cells..., compress=true, append=false)

        vtk[@Name(attribute)] = attribute;

        vtk_save(vtk)
    end

end

#function main()
#    Np = 1280
#    points = rand(2,Np);
#    filename = "polydata"
#    attribute =  [randn() for i = 1, j = 1:Np];
#    @time filenames = write_vtp(filename, points, attribute)
#end

#main()

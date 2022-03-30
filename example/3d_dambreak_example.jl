#using PostSPH

# Load all data in
cd(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreak_out\data")
rhop_array = PostSPH.readBi4Array(PostSPH.Rhop)
pos_array  = PostSPH.readBi4Array(PostSPH.Points)
vel_array  = PostSPH.readBi4Array(PostSPH.Vel)
idp_array  = PostSPH.readBi4Array(PostSPH.Idp)


Npok = PostSPH.readBi4_CurrentTotalParticles()
TypeOfParticle, NValues = PostSPH.readBi4_NumberOfParticles()
NTime = PostSPH.readBi4_Time()

# Save to vtk
function SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
    Base.Threads.@threads for i = 1:length(pos_array)
        SimData    = PostSPH.SaveVTK.SimData(
            Points = pos_array[i],
            Idp    = idp_array[i],
            Vel    = vel_array[i],
            Rhop   = rhop_array[i],
        )

        #Use 3d glyph filter in Paraview with 2d glyphs!
        PostSPH.SaveVTK.write_vtp("SimData_" * lpad(string(i), 4, "0"), SimData)
    end
end

@time SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
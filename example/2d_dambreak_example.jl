using PostSPH

# Load all data in
cd(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreakVal2D_out\data")
rhop_array = readBi4Array(PostSPH.Rhop)
pos_array = readBi4Array(PostSPH.Points)
vel_array = readBi4Array(PostSPH.Vel)
idp_array = readBi4Array(PostSPH.Idp)


Npok = PostSPH.readBi4_CurrentTotalParticles()
TypeOfParticle, NValues = PostSPH.readBi4_NumberOfParticles()
NTime = PostSPH.readBi4_Time()

# Save to vtk
for i = 1:201
    SimData101 = PostSPH.SaveVTK.SimData(
        Points = pos_array[i],
        Idp = idp_array[i],
        Vel = vel_array[i],
        Rhop = rhop_array[i],
    )

    #Use 3d glyph filter in Paraview with 2d glyphs!
    @time PostSPH.SaveVTK.write_vtp("SimData_" * lpad(string(i), 4, "0"), SimData101)
end

using PostSPH

# Load all data in
cd(raw"D:\DualSPHysics_v5.0\examples\main\03_MovingSquare\CaseMovingSquare_out\data")
rhop_array = readBi4Array(PostSPH.Rhop)
pos_array  = readBi4Array(PostSPH.Points)
vel_array  = readBi4Array(PostSPH.Vel)
idp_array  = readBi4Array(PostSPH.Idp)

# Save to vtk
SimData101 = PostSPH.SaveVTK.SimData(Points = pos_array[101],
                                     Idp    = idp_array[101],
                                     Vel    = vel_array[101],
                                     Rhop   = rhop_array[101])

#Use 3d glyph filter in Paraview with 2d glyphs!
PostSPH.SaveVTK.write_vtp("SimData101",SimData101)

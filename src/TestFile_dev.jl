#Test PostSPH in development

cd(raw"D:\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreakVal2D_Sim\particles")

#Extracts all velocity data from all .vtk files and stores in array
#Second argument is velocity, clear a few lines further down
@time data = PostSPH.readVtkArray("PartFluid",PostSPH.Cat(0))
#To see number of particles in each .vtk file do then:
@time n   = PostSPH.readVtkParticles("PartFluid")

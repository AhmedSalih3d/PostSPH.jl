#Test PostSPH in development

using BenchmarkTools

#data / particles
cd(raw"D:\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreakVal2D_Sim\data")

#Extracts all velocity data from all .vtk files and stores in array
#Second argument is velocity, clear a few lines further down
#@time data = PostSPH.readVtkArray("PartFluid",PostSPH.Cat(0))
#To see number of particles in each .vtk file do then:
#@time n   = PostSPH.readVtkParticles("PartFluid")

@time rho_arr = PostSPH.readBi4Array(PostSPH.Rhop,PostSPH._dirFiles()[1])

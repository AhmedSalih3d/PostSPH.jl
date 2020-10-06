#Test PostSPH in development

using BenchmarkTools

#data / particles
cd(raw"D:\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreakVal2D_Sim\data")

#Extracts all velocity data from all .vtk files and stores in array
#Second argument is velocity, clear a few lines further down
#@time data = PostSPH.readVtkArray("PartFluid",PostSPH.Cat(0))
#To see number of particles in each .vtk file do then:
#@time n   = PostSPH.readVtkParticles("PartFluid")

#@time rho_arr = PostSPH.readBi4Array(PostSPH.Points,PostSPH._dirFiles())
#@btime rho_arr2 = PostSPH2.readBi4Array($PostSPH2.Points,$PostSPH2._dirFiles())
#nParticles = PostSPH.readBi4Particles()


#files = PostSPH._dirFiles()[1:10]
function save_files()
    files = PostSPH._dirFiles()
    for file in files
        pos_arr = PostSPH.readBi4Array(PostSPH.Points,file)
        idp_arr = PostSPH.readBi4Array(PostSPH.Idp,file)
        vel_arr = PostSPH.readBi4Array(PostSPH.Vel,file)
        rho_arr = PostSPH.readBi4Array(PostSPH.Idp,file);
        sim_arr = PostSPH.SaveVTK.SimData(Points = pos_arr[1], Rhop=rho_arr[1], Vel=vel_arr[1])

        save_dir = "D:\\DualSPHysics_v5.0\\examples\\main\\01_DamBreak\\CaseDambreakVal2D_Sim\\particlesJulia\\"
        PostSPH.SaveVTK.write_vtp(save_dir*file[1:end-4],sim_arr)
        println("Saved ", file)
    end
end




function save_files2()
    files = PostSPH2._dirFiles()
    for file in files
        pos_arr = PostSPH2.readBi4Array(PostSPH2.Points,file)
        idp_arr = PostSPH2.readBi4Array(PostSPH2.Idp,file)
        vel_arr = PostSPH2.readBi4Array(PostSPH2.Vel,file)
        rho_arr = PostSPH2.readBi4Array(PostSPH2.Idp,file);
        sim_arr = PostSPH2.SaveVTK2.SimData(Points = pos_arr[1], Rhop=rho_arr[1], Vel=vel_arr[1])

        save_dir = "D:\\DualSPHysics_v5.0\\examples\\main\\01_DamBreak\\CaseDambreakVal2D_Sim\\particlesJulia\\"
        PostSPH2.SaveVTK2.write_vtp(save_dir*file[1:end-4],sim_arr)
    end
end

@btime save_files()
@btime save_files2()

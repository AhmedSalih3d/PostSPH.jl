#using PostSPH

# Load all data in
cd(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreak_out\data")
rhop_array = readBi4Array(PostSPH.Rhop)
pos_array  = readBi4Array(PostSPH.Points)
vel_array  = readBi4Array(PostSPH.Vel)
idp_array  = readBi4Array(PostSPH.Idp)


Npok = readBi4_CurrentTotalParticles()
TypeOfParticle, NValues = readBi4_NumberOfParticles()
NTime = readBi4_Time()

Fluid_Ids  = readBi4_Head()[4]
Column_Ids = readBi4_Head()[1]

# Save to vtk
function SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
    Base.Threads.@threads for i = 1:length(pos_array)


        Fluid_act_id    = findall(Fluid_Ids["IdRangeJulia"][1]  .<= idp_array[i] .<= Fluid_Ids["IdRangeJulia"][end])  #filter did not work?
        Column_act_id   = findall(Column_Ids["IdRangeJulia"][1] .<= idp_array[i] .<= Column_Ids["IdRangeJulia"][end]) #filter(x -> Column_Ids["IdRangeJulia"][1] <= x < Column_Ids["IdRangeJulia"][end], idp_array[i])

        SimDataFluid    = SaveVTK.SimData(
            Points = pos_array[i],
            Idp    = idp_array[i],
            Vel    = vel_array[i],
            Rhop   = rhop_array[i],
        )

        SimDataColumn   = SaveVTK.SimData(
            Points = pos_array[i],
            Idp    = idp_array[i],
            Vel    = vel_array[i],
            Rhop   = rhop_array[i],
        )

        #Use 3d glyph filter in Paraview with glyphs!
        SaveVTK.write_vtp("SimData_FLUID_" * lpad(string(i), 4, "0"), SimDataFluid,Fluid_act_id)
        SaveVTK.write_vtp("SimData_COLUMN_" * lpad(string(i), 4, "0"), SimDataColumn,Column_act_id)
    end
end

@time SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
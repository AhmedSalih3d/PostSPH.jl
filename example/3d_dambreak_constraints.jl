using PostSPH

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
Column_Ids = readBi4_Head()[2]

# Save to vtk
function SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
    Base.Threads.@threads for i = 1:length(pos_array)


        Fluid_act_id    = findall(Fluid_Ids["Begin"] .<= idp_array[i] .<=  Fluid_Ids["Begin"] +  Fluid_Ids["Count"]-1)  
        Column_act_id   = findall(Column_Ids["Begin"] .<= idp_array[i] .<= Column_Ids["Begin"]+ Column_Ids["Count"]-1)

        SimDataFluid    = PostSPH.SaveVTK.SimData(
            Points = pos_array[i],
            Idp    = idp_array[i],
            Vel    = vel_array[i],
            Rhop   = rhop_array[i],
        )

        SimDataColumn   = PostSPH.SaveVTK.SimData(
            Points = pos_array[i],
            Idp    = idp_array[i],
            Vel    = vel_array[i],
            Rhop   = rhop_array[i],
        )

        # Bit difficult to comprehend, but it broadcasts a constraint function, on a parameter (constraint array) to determine how to
        # constrain a selected array
        function Constrain(ArrayToConstrain::AbstractArray,ConstraintArray::AbstractArray,ConstraintFunction::Function)
            return ConstraintFunction.(ConstraintArray[ArrayToConstrain])
        end

        fx(x) = 0.2 <= x <= 1.3
        fy(x) = 0.1 <= x <= 0.5
        fxc   = Constrain(Fluid_act_id,pos_array[i][1:3:end],fx)
        fyc   = Constrain(Fluid_act_id,pos_array[i][2:3:end],fy)

        fc_id = fxc .* fyc

        #Use 3d glyph filter in Paraview with glyphs!
        PostSPH.SaveVTK.write_vtp(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreak_out\particles_vtp\SimData_FLUID_"  * lpad(string(i-1), 4, "0"), SimDataFluid,Fluid_act_id[fc_id]) #[fc_id]
        PostSPH.SaveVTK.write_vtp(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreak_out\particles_vtp\SimData_COLUMN_" * lpad(string(i-1), 4, "0"), SimDataColumn,Column_act_id)
    end
end

@time SaveVTK_out(pos_array,idp_array,vel_array,rhop_array)
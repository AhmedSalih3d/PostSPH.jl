using Plots

plot(0)
title!("SumKE")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\GenericWaveTank_mDBC.xml")
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
res = map(x->0.5*sum(sqrt.(sum(x.^2,dims=1)).^2),p2)

plot!(res,label="Z1_SumKE")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z3\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z3\GenericWaveTank_mDBC.xml")
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
res = map(x->0.5*sum(sqrt.(sum(x.^2,dims=1)).^2),p2)

plot!(res,label="Z3_SumKE")


cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\03_NoPlate\GenericWaveTank_mDBC_NoPlate1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\03_NoPlate\GenericWaveTank_mDBC_NoPlate1\GenericWaveTank_mDBC_postedit.xml")
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
res = map(x->0.5*sum(sqrt.(sum(x.^2,dims=1)).^2),p2)
plot!(res,label="Base_SumKe"


cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\data")

bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
res = map(x->0.5*sum(sqrt.(sum(x.^2,dims=1)).^2),p2)
plot!(res,label="Z2")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\GenericWaveTank_mDBC_out\data")

bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\GenericWaveTank_mDBC_out\GenericWaveTank_mDBC.xml")
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
res = map(x->0.5*sum(sqrt.(sum(x.^2,dims=1)).^2),p2)
plot!(res,label="Z2L")

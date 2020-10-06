using Plots

plot(0)
title!("SumKE - For particles in 0 to 1.5 meters in x")
xlabel!("Moving Average Sample - 100")
ylabel!("SumKE")

na = 100
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\GenericWaveTank_mDBC.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res = zeros(length(p2))
for i in eachindex(p2)
    res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
end

plot!(moving_average(res,na),label="Z1")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\data")

bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res = zeros(length(p2))
for i in eachindex(p2)
    res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
end

plot!(moving_average(res,na),label="Z2")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z3\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z3\GenericWaveTank_mDBC.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res = zeros(length(p2))
for i in eachindex(p2)
    res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
end

plot!(moving_average(res,na),label="Z3")


cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\03_NoPlate\GenericWaveTank_mDBC_NoPlate1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\03_NoPlate\GenericWaveTank_mDBC_NoPlate1\GenericWaveTank_mDBC_postedit.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res = zeros(length(p2))
for i in eachindex(p2)
    res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
end

plot!(moving_average(res,na),label="NoPlate")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\05_L\GenericWaveTank_mDBC_Z2L\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\05_L\GenericWaveTank_mDBC_Z2L\GenericWaveTank_mDBC.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res = zeros(length(p2))
for i in eachindex(p2)
    res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
end

plot!(moving_average(res,na),label="Z2L",legend=:outertopright)

# cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\GenericWaveTank_mDBC_out\data")
# bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\GenericWaveTank_mDBC_out\GenericWaveTank_mDBC.xml")
# p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
# indices =  map(x-> 0 .<= x[1,:] .<= 1.5,p1)
# p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
#
# res = zeros(length(p2))
# for i in eachindex(p2)
#     res[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices[i]].^2,dims=1)).^2)
# end
#
# plot!(moving_average(res,na),label="OUT",legend=:outertopright)

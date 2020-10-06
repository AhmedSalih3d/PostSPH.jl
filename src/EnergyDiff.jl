using Plots

plot(0)
title!("SumKE - Difference between \n left and right edge for particles in \n -0.2 to 0 and 1.5 to 1.7 meters in x")
xlabel!("Moving Average Sample - 200")
ylabel!("SumKE")

na = 100
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z1\GenericWaveTank_mDBC.xml")
p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices1 =  map(x-> -0.2 .<= x[1,:] .<= 0,p1)
indices2 =  map(x-> 1.5 .<= x[1,:] .<= 1.7,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res1 = zeros(length(p2))
res2 = zeros(length(p2))
for i in eachindex(p2)
    res1[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices1[i]].^2,dims=1)).^2)
    res2[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices2[i]].^2,dims=1)).^2)
end

plot!(moving_average(res1,na) .- moving_average(res2,na),label="Z1_Left-Right")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")

p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices1 =  map(x-> -0.2 .<= x[1,:] .<= 0,p1)
indices2 =  map(x-> 1.5 .<= x[1,:] .<= 1.7,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res1 = zeros(length(p2))
res2 = zeros(length(p2))
for i in eachindex(p2)
    res1[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices1[i]].^2,dims=1)).^2)
    res2[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices2[i]].^2,dims=1)).^2)
end

plot!(moving_average(res1,na) .- moving_average(res2,na),label="Z2_Left-Right")


cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z3\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")

p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices1 =  map(x-> -0.2 .<= x[1,:] .<= 0,p1)
indices2 =  map(x-> 1.5 .<= x[1,:] .<= 1.7,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res1 = zeros(length(p2))
res2 = zeros(length(p2))
for i in eachindex(p2)
    res1[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices1[i]].^2,dims=1)).^2)
    res2[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices2[i]].^2,dims=1)).^2)
end

plot!(moving_average(res1,na) .- moving_average(res2,na),label="Z3_Left-Right")

cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\03_NoPlate\GenericWaveTank_mDBC_NoPlate1\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")

p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices1 =  map(x-> -0.2 .<= x[1,:] .<= 0,p1)
indices2 =  map(x-> 1.5 .<= x[1,:] .<= 1.7,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res1 = zeros(length(p2))
res2 = zeros(length(p2))
for i in eachindex(p2)
    res1[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices1[i]].^2,dims=1)).^2)
    res2[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices2[i]].^2,dims=1)).^2)
end

plot!(moving_average(res1,na) .- moving_average(res2,na),label="NoPlate1_Left-Right")


cd(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\05_L\GenericWaveTank_mDBC_Z2L\data")
bodies = MkArray(raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\02_ChangingPlateZ\GenericWaveTank_mDBC_Z2\GenericWaveTank_mDBC.xml")

p1 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
indices1 =  map(x-> -0.2 .<= x[1,:] .<= 0,p1)
indices2 =  map(x-> 1.5 .<= x[1,:] .<= 1.7,p1)
p2 = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)

res1 = zeros(length(p2))
res2 = zeros(length(p2))
for i in eachindex(p2)
    res1[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices1[i]].^2,dims=1)).^2)
    res2[i] = 0.5*sum(sqrt.(sum(p2[i][:,indices2[i]].^2,dims=1)).^2)
end

plot!(moving_average(res1,na) .- moving_average(res2,na),label="Z2L_Left-Right",legend=:outertopright)

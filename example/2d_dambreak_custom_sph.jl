using PostSPH
using CellListMap
using StaticArrays

##
const dp = 0.01
const rhop0 = 1000
const H = 0.014142135624
const Hinv = 1/H
const cutoff = 2*H
const mass = rhop0*dp^2
const c0   = 88.36718777351693
const b    = 1115537.1429
const gamma= 7
const aD   = 7/(4*π*H^2)
const alpha = 0.1

const G = SVector{2,Float64}(0,-9.81)

dt = 0.00001

function EquationOfState(rho,rhop0,b,gamma)
    return b*((rho/rhop0)^gamma - 1)
end

function WendlandDerivative(q,aD)
    if q > 2
        val = 0
    else
        val = aD*((5.0/8.0)*q*(q-2.0)^3)
    end
    return val
end

# Load all data in
cd(raw"C:\Users\Ahmed Salih\Documents\DualSPHysics_v5.0\examples\main\01_DamBreak\CaseDambreakVal2D_out\data")
rhop_array = readBi4Array(PostSPH.Rhop,"Part_0000.bi4")[1]
pos_array  = readBi4Array(PostSPH.Points,"Part_0000.bi4")[1]
vel_array  = readBi4Array(PostSPH.Vel,"Part_0000.bi4")[1]
idp_array  = readBi4Array(PostSPH.Idp,"Part_0000.bi4")[1]

const N    = length(idp_array)

Wall_Ids   = readBi4_Head()[1]
Fluid_Ids  = readBi4_Head()[2]
Fluid_act_Ids   = findall(Fluid_Ids["Begin"] .<= idp_array .<=Fluid_Ids["Begin"]+Fluid_Ids["Count"]-1)
Wall_act_Ids   = findall(Wall_Ids["Begin"] .<= idp_array .<= Wall_Ids["Begin"]+ Wall_Ids["Count"]-1)

ps_SA = Vector{SVector{2,Float64}}(undef,N);
x_view = @view pos_array[1:3:end]
y_view = @view pos_array[2:3:end]
z_view = @view pos_array[3:3:end]
for i = 1:N
    x__      = x_view[i]
    z__      = z_view[i]
    ps_SA[i] = SVector{2,Float64}(x__,z__)
end

vs_SA = Vector{SVector{2,Float64}}(undef,N);
vx_view = @view vel_array[1:3:end]
vy_view = @view vel_array[2:3:end]
vz_view = @view vel_array[3:3:end]
for i = 1:N
    vx__      = vx_view[i]
    vz__      = vz_view[1]
    vs_SA[i] = SVector{2,Float64}(vx__,vz__)
end

#preallocate
pressure = [ 0.0 for i in 1:N ]
rhop     = Float64.(deepcopy(rhop_array))
rhop_intermediate = deepcopy(rhop)
pos      = deepcopy(ps_SA)
vel      = deepcopy(vs_SA)
vel0     = deepcopy(vs_SA) .* 0
vel_intermediate = deepcopy(vel)


function calc_continuity_equation!(x, y, i, j, d2, Hinv,aD,mass,rhop,vel) 
    d  = sqrt(d2)+1e-6 #To avoid division by zero
    q  = d*Hinv

    Wq = WendlandDerivative(q,aD)

    Wg = Wq * ((x-y)/d) * Hinv

    drhoi = (mass/rhop[j]) * sum(vel[i] .* Wg)

    rhop[i] += drhoi
    rhop[j] -= drhoi

    return rhop
end

function calc_momentum_equation!(x, y, i, j, d2, H,Hinv,alpha,c0,aD,mass,rhop,pressure,vel,G) 
    d  = sqrt(d2)+1e-6
    q  = d*Hinv

    Wq = WendlandDerivative(q,aD)

    Wg = Wq * ((x-y)/d) * Hinv

    Pa = pressure[i]
    Pb = pressure[j]
    rhoab = rhop[i]*rhop[j]

    # velpos = (vel[i]-vel[j]) .* (x-y)
    # cond = sum(velpos)
    # if cond < 0
    #     μ = H * cond / (d + 0.001*H^2)  #Weird in documentation!
    #     Π =  (-alpha*c0*μ)/((rhop[i]+rhop[j])/2)
    # else
    #     Π = 0
    # end

    dvi = -mass*(((Pa+Pb)/rhoab))*Wg+G;

    vel[i] = vel[i]+dvi
    vel[j] = vel[j]-dvi


    return vel
end

# Run pairwise computation
# map_pairwise!(
#     (x, y, i, j, d2, pressure) -> calc_pressure!(x, y, i, j, d2, b, rhop0, gamma, rhop,pressure),
#     pressure,box,cl,parallel=false
# )

cd(raw"C:\Users\Ahmed Salih\Documents\RunSim")

for time_step = 1:length(range(0,2,step=dt))

    box = Box(limits(pos),cutoff)
    cl  = CellList(pos,box)

    pressure = EquationOfState.(rhop,rhop0,b,gamma)

    map_pairwise!(
        (x, y, i, j, d2, rhop_intermediate) -> calc_continuity_equation!(x, y, i, j, d2, Hinv,aD,mass,rhop_intermediate,vel_intermediate),
        rhop_intermediate,box,cl,parallel=false
    )

    map_pairwise!(
        (x, y, i, j, d2, vel_intermediate) -> calc_momentum_equation!(x, y, i, j, d2, H,Hinv,alpha,c0,aD,mass,rhop,pressure,vel0,G),
        vel_intermediate,box,cl,parallel=false
    )

    # To ensure wall does not move

    vel[Wall_act_Ids] *= 0;
    ##

    if mod(time_step,1) == 0
    NN = length(pos)
    fp = first.(pos)
    lp = last.(pos)
    fv = first.(vel)
    lv = last.(vel)
    pos_vtk = zeros(Float64,NN*3)
    vel_vtk = zeros(Float64,NN*3)
    for ii = 1:NN-1
        pos_vtk[3*(ii-1)+1]   = fp[ii]
        pos_vtk[3*(ii-1)+2] = 0
        pos_vtk[3*(ii-1)+3] = lp[ii]

        vel_vtk[3*(ii-1)+1]   = fv[ii]
        vel_vtk[3*(ii-1)+2] = 0
        vel_vtk[3*(ii-1)+3] = lv[ii]
    end

    SimData    = PostSPH.SaveVTK.SimData(
                Points = pos_vtk,
                Idp    = idp_array,
                Vel    = vel_vtk,
                Rhop   = rhop,
            )

    PostSPH.SaveVTK.write_vtp("SimData_Fluid" * lpad(string(time_step-1), 4, "0"), SimData,Fluid_act_Ids)
    PostSPH.SaveVTK.write_vtp("SimData_Wall" * lpad(string(time_step-1), 4, "0"), SimData,Wall_act_Ids)
    end
    println(time_step)

    # Time Stepping
    rhop = rhop + 2*dt*rhop_intermediate
    vel  = vel  + 2*dt*vel_intermediate
    pos  = pos + dt*vel + 0.5*dt^2*vel_intermediate

    rhop_intermediate = rhop
    vel_intermediate = vel
end



# # Save to vtk
# for i = 1:201
#     SimData    = PostSPH.SaveVTK.SimData(
#         Points = pos_array[i],
#         Idp    = idp_array[i],
#         Vel    = vel_array[i],
#         Rhop   = rhop_array[i],
#     )

#     #Use 3d glyph filter in Paraview with 2d glyphs!
#     @time PostSPH.SaveVTK.write_vtp("SimData_" * lpad(string(i), 4, "0"), SimData)
# end

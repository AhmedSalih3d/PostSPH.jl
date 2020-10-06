include("PostSPH.jl")
include("PistonForce.jl")
using Statistics #to use "mean"
using Plots
using PlotThemes
gr()


theme(:vibrant) #vibrant
default(linewidth = 1.5) #default() to reset
function mag(arr)
    return sqrt.(arr[1,:].^2 .+ arr[3,:].^2)
end


#top_path  = raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\08_DpSPS\GenericWaveTank_mDBC_1.2_Dp2"
top_path  = raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\09_DpART\03_GenericWaveTank_mDBC_1.2_Dp1ART"
mass_bound = 0.4 #0.4
ke1,pe1,ke2,pe2,de1,de2 = energy_calc(top_path,mass_bound)
gr_plot_res(top_path,ke1,pe1,ke2,pe2,de1,de2)

function energy_calc(top_path,mass_bound)
    data_path = raw"data"
    name_xml  = raw"GenericWaveTank_mDBC.xml"
    cd(joinpath(top_path,data_path))
    bodies = MkArray(joinpath(top_path,name_xml))

    ## FLuid

    fluid_pos = PostSPH.readBi4Body(bodies[4][1],PostSPH.Points)
    fluid_vel = PostSPH.readBi4Body(bodies[4][1],PostSPH.Vel)
    g = 9.80675
    ir = 6 #influence range
    pl = 1.5 #plate length
    orig = 0 #start point of plate
    d  = 1.2 #depth
    indices1 =  map(x-> -ir .<= x[1,:] .<= orig,fluid_pos)
    indices2 =  map(x-> pl .<= x[1,:] .<= pl + ir,fluid_pos)

    function ke_all(fluid_vel,indices1,indices2)
        ke1 = zeros(size(fluid_vel)[1])
        ke2 = zeros(size(fluid_vel)[1])
        pe1 = zeros(size(fluid_vel)[1])
        pe2 = zeros(size(fluid_vel)[1])


        for i in eachindex(fluid_vel)

            ke = 0.5 * mass_bound * sum(mag(fluid_vel[i][:,indices1[i]]).^2)

            pe = mass_bound .* g * sum((fluid_pos[i][3,indices1[i]] .+ d))

            ke1[i] = ke
            pe1[i] = pe

        end

        for i in eachindex(fluid_vel)

            ke = 0.5 * mass_bound * sum(mag(fluid_vel[i][:,indices2[i]]).^2)

            pe = mass_bound * g * sum((fluid_pos[i][3,indices2[i]] .+ d))

            ke2[i] = ke
            pe2[i] = pe

        end

        return ke1,pe1,ke2,pe2
    end

    ke1,pe1,ke2,pe2 = ke_all(fluid_vel,indices1,indices2)

    function de(Work, te1,te2)
        n = length(Work)

        te1_ini = te1[1]
        te2_ini = te2[1]

        de1 = zeros(n)
        de2 = zeros(n)
        for i in 1:n-1
            #de1[i+1] = Work[i+1] + te1_ini - te1[i+1]
            #de2[i+1] = Work[i+1] + te2_ini - te2[i+1]
            #de1[i+1] = Work[i+1] + te1[i] - te1[i+1]
            #de2[i+1] = Work[i+1] + te2[i] - te2[i+1]
            de1[i+1] = Work[i] + te1[i] - te1[i+1]
            de2[i+1] = Work[i] + te2[i] - te2[i+1]
        end

        return de1,de2
    end

    te1 = ke1 .+ pe1
    te2 = ke2 .+ pe2
    de1,de2 = de(Work,te1,te2)

     #df = DataFrame([t,ke1,pe1,ke2,pe2,de1,de2],[:t,:ke1,:pe1,:ke2,:pe2,:de1,:de2])
     #CSV.write("path",df)
    return ke1,pe1,ke2,pe2,de1,de2
end

function gr_plot_res(top_path,ke1,pe1,ke2,pe2,de1,de2)
    #Plots

    case_name = basename(top_path)

    p0  = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(t,F,label="Piston Force")
    xlabel!("Time [s]")
    ylabel!("Force [N]")
    title!("2D Force exerted by Piston\n$case_name")

    p00 =  plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(t,Work,label="Piston Work")
    xlabel!("Time [s]")
    ylabel!("Energy [W]")
    title!("Work Produced by Piston\n$case_name")

    p1 = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(t,ke1,label="Kinetic Energy Left of Plate")
    plot!(t,ke2,label="Kinetic Energy Right of Plate")
    xlabel!("Time [s]")
    ylabel!("Energy [W]")
    title!("Kinetic Energy\n$case_name")

    p2 = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(t,pe1,label="Potential Energy Left of Plate")
    plot!(t,pe2,label="Potential Energy Right of Plate")
    xlabel!("Time [s]")
    ylabel!("Energy [W]")
    title!("Potential Energy\n$case_name")

    p3 = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
    plot!(t,de1,label="Dissipation of Energy Left of Plate")
    plot!(t,de2,label="Dissipation of Energy Right of Plate")
    xlabel!("Time [s]")
    ylabel!("Energy [W]")
    title!("Dissipation Energy\n$case_name")

    # Display plots

    pdf_path  = mkpath(joinpath(top_path,"pdf"))

    display(p0)
    savefig(p0, joinpath(pdf_path,"gr_FORCE_$case_name.pdf"))
    sleep(1)
    display(p00)
    savefig(p00, joinpath(pdf_path,"gr_WORK_$case_name.pdf"))
    sleep(1)
    display(p1)
    savefig(p1, joinpath(pdf_path,"gr_KE_$case_name.pdf"))
    sleep(1)
    display(p2)
    savefig(p2, joinpath(pdf_path,"gr_PE_$case_name.pdf"))
    sleep(1)
    display(p3)
    savefig(p3, joinpath(pdf_path,"gr_DE_$case_name.pdf"))
    sleep(1)

    return "Plotting done"
end

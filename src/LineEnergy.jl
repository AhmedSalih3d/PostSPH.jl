using Statistics #to use "mean"
using Plots
using DataFrames
using CSV
using Gtk
using PGFPlotsX
using Colors
pgfplotsx()

# Pop redundant entries
for _ in eachindex(PGFPlotsX.CUSTOM_PREAMBLE)
    pop!(PGFPlotsX.CUSTOM_PREAMBLE)
end
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{pgfplotstable}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{booktabs,colortbl}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{siunitx}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage{amsmath}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usetikzlibrary{shapes.geometric}")
push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepackage[obeyspaces]{url}")


## Extended energy calc

function read_pres(top_path)
    pres_file = raw"EF10pts_01_GenericWaveTank_mDBC_1_Press.csv"
    fname = joinpath(joinpath(top_path,"energyflux"),pres_file)
    file = open(fname,"r")

    f_parse_f64(f) = parse.(Float64,split(f,";")[3:end])

    f_posx   = readline(file)
    f_posy   = readline(file)
    f_posz   = readline(file)
    f_header = readline(file)

    P_posx   = f_parse_f64(f_posx)
    P_posy   = f_parse_f64(f_posy)
    P_posz   = f_parse_f64(f_posz)
    P_header = replace.(split(f_header,";"),Ref(r" \[m\]"=>""))
    P_header = replace.(P_header,Ref(r" \[s\]"=>""))



    P_df = CSV.read(fname,DataFrame; delim=";",skipto=5,header=false)
    rename!(P_df,Symbol.(P_header))

    close(file)
    return P_posx,P_posz,P_header,P_df
end

function read_vel(top_path)
    vel_file  = raw"EF10pts_01_GenericWaveTank_mDBC_1_Vel.csv"
    fname = joinpath(joinpath(top_path,"energyflux"),vel_file)
    file = open(fname,"r")

    f_parse_f64(f) = parse.(Float64,split(f,";")[3:end])

    f_pos    = readline(file)
    f_header = readline(file)

    f_pos_parsed = f_parse_f64(f_pos)

    V_pos = reshape(f_pos_parsed,(3,Int(length(f_pos_parsed)/3)))

    V_header = replace.(split(f_header,";"),Ref(r" \[m\]"=>""))
    V_header = replace.(V_header,Ref(r" \[s\]"=>""))
    V_header = replace.(V_header,Ref(r" \[m/s\]"=>""))

    V_header_x = findall( x -> occursin("x", x), V_header)

    V_df = CSV.read(fname,DataFrame; delim=";",skipto=3,header=false)
    rename!(V_df,Symbol.(V_header))

    V_df_x = V_df[!,V_header_x]

    close(file)
    return V_pos,V_header_x,V_df_x
end

# Due to only using x-velocities the ids for pressure also work for velocity Alhamdu lillah
function unique_ids(P_posx)
    u_id = unique(P_posx)

    arr_id = fill(BitArray{1}(),length(u_id))
    for i = 1:length(u_id)
        arr_id[i] = P_posx .== u_id[i]
    end
    return arr_id
end

function calc_energy(P_posz,P_header,P_df,V_df_x,ids,no)
    # Split pressure dataframe in part / time and pressure values - therefore also shrink header
    P_df_PT  = P_df[!,P_header[1:2]]
    P_df     = select(P_df, Not(P_header[1:2]))
    P_header = P_header[3:end] 
    # Define ω -  hard-coded
    T = 2
    ω = 2pi/T
    # Define gravity and density - hard-coded
    g   = 9.80675
    rho = 1000
    H   = -1.2
    # Preallocate e_final
    e_final = fill(Array{Float64,1}(),size(P_df)[1])
    #e_final = zeros(size(P_df)[1])
    # Extract gauge pressure;
    p_g = P_df[P_header[ids[no]]]
    # Extract velocity in x
    V_x  = V_df_x[ids[no]]
    # Multiply p_d with relevant x_velocities
    #  P_posz == V_pos[3,:] -> true
    # -(pg_r*H .- rho*g*posz.^2 / 2) -> always positive
    posz = P_posz[ids[no]]
    @inbounds for i = 1:length(e_final)
        pg_r = convert(Array,p_g[i,:])
        Vx_r = convert(Array,V_x[i,:])
        e_z = -(pg_r*H .- rho*g*posz.^2 / 2) .* Vx_r
        e_t = 2pi/ω * e_z
        e_final[i] = ω/2pi * e_t
        #e_final[i] = ω/2pi * sum(e_t)
    end

    return posz,e_final,V_x
end

folder = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER,select_multiple=false)

P_posx,P_posz,P_header,P_df = read_pres(folder)
V_pos,V_header_x,V_df_x     = read_vel(folder)
ids = unique_ids(P_posx)

p = plot()
for i = 15
    posz,e_final,V_x = calc_energy(P_posz,P_header,P_df,V_df_x,ids,i)
    plot!(posz,e_final)
    #println(V_x);
end
display(p)
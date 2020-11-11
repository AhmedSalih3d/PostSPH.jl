using Statistics #to use "mean"
using Plots
using Plots.Measures
using DataFrames
using CSV
using Gtk
using PGFPlotsX
using Colors
using LaTeXStrings
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


theme(:vibrant,
             fglegend = colorant"#000000", #legend
             guidefontcolor = colorant"#000000",
             showaxis = true,
             #formatter = :scientific,
             aspect_ratio = :equal, #:equal
             #widen=true,
             #ylims=(0.1,100),xlims=(0,100),
             thickness_scaling = 1, 
             framestyle=:box,
             size=(1300,400),
             xtickfontsize=22,
             ytickfontsize=22,
             yguidefontsize=22,
             xguidefontsize=22,
             legendfontsize=22,
             titlefontsize=22,
             linewidth = 1,
             margin=10mm,
               )

function p64(string)
    parse(Float64,string)
end

macro name(arg)
   string(arg)
end

# function to read specific line
function readline_n(file,n)
        open(file, "r") do io
            for i = 1:n
                if i == n
                    return readline(io)
                    break
                end
                readline(io)
            end
        end
    end

# Function to load all exported case variables

function extract_vars_from_line!(dict,d)
    dd = split(d,":")[2]
    ddd = split(dd," ")
    indices = findall(map(!,isempty.(ddd)))
    dddd = ddd[indices]
    n    = length(dddd)
    ddddd = split.(dddd,"=")

    for i = 1:n
        d   = ddddd[i]
        s   = d[1]
        v   = strip(d[2],['[',']'])
        push!(dict, Symbol(s) => v)
    end
end

function case_vars(path)
    fpath = joinpath(path,"Run.out")
    # User vars + cte
    d =  readline_n(fpath,68)
    # Some parameters
    e =  readline_n(fpath,69)

    dict = Dict{Symbol,String}()

    extract_vars_from_line!(dict,d)
    extract_vars_from_line!(dict,e)

    return dict

end


searchdir(path,key) = filter(x->occursin(key,x), readdir(path))


function plot_it(path,fname,posx,col_mean,mht)
    function table_tex(title_jl_to_use,posx)
        #title_jl_to_use = replace(path, "_" => "\\_")
        "
        \\node (title) at (\$(group c1r1.center)!0.5!(group c1r1.center)+(0,4.5 cm)\$) {\\begin{tabular}{c} \\textbf{Wave Height measured around SWL} \\\\ \\path{$title_jl_to_use} \\end{tabular}};
        %(group c1r1.south east)
        \\node [below right, xshift=-1.7mm,yshift=0mm] at (leg1.south west) {
        \\renewcommand{\\arraystretch}{1} % Default value: 1
        %\\LARGE
        \\resizebox{0.225\\textwidth}{!}{%
        \\begin{tabular}{@{}llr@{}}
        \\toprule
        \\textbf{Elevation} & \\textbf{X} \\\\ \\midrule
        0                  &  $(posx[1])                      \\\\
        1                  &  $(posx[2])                       \\\\
        2                  &  $(posx[3])                       \\\\
        3                  &  $(posx[4])                       \\\\
        4                  &  $(posx[5])                       \\\\
        5                  &  \\phantom{-}$(posx[6])                       \\\\ \\bottomrule
        \\end{tabular}%
        }
        };
        \\node [below right, xshift=-1.7mm,yshift=0mm] at (leg2.south west) {
        \\renewcommand{\\arraystretch}{1} % Default value: 1
        %\\LARGE
        \\resizebox{0.225\\textwidth}{!}{%
        \\begin{tabular}{@{}llr@{}}
        \\toprule
        \\textbf{Elevation} & \\textbf{X} \\\\ \\midrule
        6                  &  $(posx[7])                      \\\\
        7                  &  $(posx[8])                       \\\\
        8                  &  $(posx[9])                       \\\\
        9                  &  $(posx[10])                       \\\\
        10                 &  $(posx[11])                       \\\\
        11                 &  $(posx[12])                       \\\\ \\bottomrule
        \\end{tabular}%
        }
        };
        % Draw horizontal line
        \\node[below right, xshift=-50mm, yshift=12.5mm] (a) at (group c1r3.north west) {};
        \\node[below right, xshift=50mm, yshift=12.5mm] (b) at (group c1r3.north east) {};
        \\draw[line width=0.625mm,   black] (a)--(b);
        "
    end

    tp = @pgf TikzPicture()

    hlines     = @pgf [HLine({color = "gray", style ="{dashed}"},0.3/2), HLine({color = "gray", style ="{dashed}"}, -0.3/2)]

    cc = distinguishable_colors(7)

    cc = [
        colorant"rgb(0,0,0)"
        colorant"rgb(135,135,135)";
        colorant"rgb(238,127,0)";
        colorant"rgb(139,173,63)";
        colorant"rgb(226,0,122)";
        colorant"rgb(0,171,164)";
        colorant"rgb(0,61,115)";
        ]

    p1 = @pgf GroupPlot(
    {
    group_style =
    {
    group_size="1 by 2",
    #yticklabels_at="edge left",
    #xticklabels_at="edge bottom",
    vertical_sep="45pt",
    },
    no_markers
    },
    #Height Mean
    {legend_pos  = "outer north east",title=raw"Height Mean - $\langle h \rangle_t$", xlabel = "Elevation X", ylabel = raw"Height Mean - $\left[\SI{}{\meter}\right]$", height = "8.25cm", width = "15cm",grid,legend_style="{name=leg5,draw=none}",xtick=posx[1:25:end],yshift=-25},
    PlotInc({color = cc[1],"ultra thick"},Table(posx, col_mean)),
    LegendEntry(raw"Height Mean"),
    #Height variance
    {legend_pos  = "outer north east",title=raw"Height Variance - $\langle \left( h - \langle h \rangle \right)^2 \rangle_t$", xlabel = "Elevation X", ylabel = raw"Height Variance - $\left[\SI{}{\meter\squared}\right]$", height = "8.25cm", width = "15cm",grid,legend_style="{name=leg6,draw=none}",xtick=posx[1:25:end],yshift=-25},
    PlotInc({color = cc[1],"ultra thick"},Table(posx, mht)),
    LegendEntry(raw"Height Variance"),
    )

    push!(tp,p1)

    pgfsave(joinpath(path,fname)*"_height"*".pdf", tp; include_preamble = true, dpi = 600)

    println("Done with $fname")
end


# function make_plots(top_path)
#     fpath = "elevation"


#     final_path = joinpath.(top_path,fpath)

#     #file_ext = "csv"
#     file_ext = "grid"
#     d = []
#     for path in final_path

#         casevar = case_vars(dirname(path))

#         csv_file = searchdir(path,file_ext)

#         true_path = joinpath(path,csv_file[1])

#         file = open(true_path,"r")

#         f_parse_f64(f) = parse.(Float64,split(f,";")[3:end])

#         f_posx = readline(file)
#         f_posy = readline(file)
#         f_posz = readline(file)
#         f_header = readline(file)

#         posx   = f_parse_f64(f_posx)
#         posy   = f_parse_f64(f_posy)
#         posz   = f_parse_f64(f_posz)
#         header = replace.(split(f_header,";"),Ref(r" \[m\]"=>""))
#         header = replace.(header,Ref(r" \[s\]"=>""))
#         n_ign  = 2
#         n_start= n_ign + 1
#         n_col  = length(header) - n_ign


#         df = CSV.read(true_path,DataFrame; delim=";",skipto=5,header=false)
#         rename!(df,Symbol.(header))

#         col_mean = mean.(eachcol(df[:,n_start:end]))

#         mht = similar(col_mean)
#         for i = 1:length(col_mean)
#             mht[i] = mean((df[!,n_ign+i] .- col_mean[i]).^2)
#         end


#         plot_it(path,csv_file[1],posx,col_mean,mht)

#         df_res = DataFrame([posx,col_mean,mht],[:posx,:col_mean,:mht])
#         CSV.write(joinpath(path,(@name df_res)*"_elevation.csv"),df_res)

#         d = df
#     end

#     return d
# end

function make_plots(top_path)
    fpath = "elevation"


    final_path = joinpath.(top_path,fpath)

    #file_ext = "csv"
    file_ext = "grid"
    d = []
    for path in final_path

        casevar = case_vars(dirname(path))

        csv_file = searchdir(path,file_ext)

        true_path = joinpath(path,csv_file[1])

        file = open(true_path,"r")

        f_parse_f64(f) = parse.(Float64,split(f,";")[3:end])

        f_posx = readline(file)
        f_posy = readline(file)
        f_posz = readline(file)
        f_header = readline(file)

        posx   = f_parse_f64(f_posx)
        posy   = f_parse_f64(f_posy)
        posz   = f_parse_f64(f_posz)
        header = replace.(split(f_header,";"),Ref(r" \[m\]"=>""))
        header = replace.(header,Ref(r" \[s\]"=>""))
        n_ign  = 2
        n_start= n_ign + 1
        n_col  = length(header) - n_ign


        df = CSV.read(true_path,DataFrame; delim=";",skipto=5,header=false)
        rename!(df,Symbol.(header))

        col_mean = mean.(eachcol(df[:,n_start:end]))

        mht = similar(col_mean)
        for i = 1:length(col_mean)
            mht[i] = mean((df[!,n_ign+i] .- col_mean[i]).^2)
        end


        plot_it(path,csv_file[1],posx,col_mean,mht)

        df_res = DataFrame([posx,col_mean,mht],[:posx,:col_mean,:mht])
        CSV.write(joinpath(path,"df_res_elevation.csv"),df_res)

        d = df
    end

    return d
end

function plot_df_res(folders,filename)
    p1 = plot(legend=:outertopright,palette=:tab20)
    for folder in folders
        df = CSV.read(joinpath(joinpath(folder,"elevation"),"df_res_elevation.csv"),DataFrame; delim=",")
        plot!(df[!,"posx"],df[!,"col_mean"],label=basename(folder))
        ylims!(-0.05,0.05)
        ylabel!("[m]")
        title!(L"Height Mean - $\langle h \rangle_t$")
    end

    p2 = plot(legend=:outertopright,palette=:tab20)
    for folder in folders
        df = CSV.read(joinpath(joinpath(folder,"elevation"),"df_res_elevation.csv"),DataFrame; delim=",")
        plot!(df[!,"posx"],df[!,"mht"],label=basename(folder))
        ylims!(0,0.03)
        xlabel!("Horizontal Position")
        ylabel!(L"[$\textrm{m}^2$]")
        title!(L"Height Variance - $\langle \left( h - \langle h \rangle \right)^2 \rangle_t$")
    end

    l = @layout grid(2,1)

    pall = plot(p1,p2, layout = l,size=(800,800))

    display(pall)

    save_path = filename*".pdf"
    savefig(save_path)
end

## Extended energy calc

function read_pres(fname)
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

function read_vel(fname)
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

function calc_energy(P_posx,P_posz,P_header,P_df,V_df_x,ids)
    # Get relevant ids to get correct vertical line - for both pressure and velocity
    ids = unique_ids(P_posx)
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
    # Extract gauge pressure;
    p_g = P_df[P_header[ids[1]]]
    # Extract velocity in x
    V_x  = V_df_x[ids[1]]
    # Multiply p_d with relevant x_velocities
    #  P_posz == V_pos[3,:] -> true
    # -(pg_r*H .- rho*g*posz.^2 / 2) -> always positive
    @inbounds for i = 1:length(e_final)
        pg_r = convert(Array,p_g[i,:])
        Vx_r = convert(Array,V_x[i,:])
        posz = P_posz[ids[1]]
        e_z = -(pg_r*H .- rho*g*posz.^2 / 2) .* Vx_r
        e_t = 2pi/ω * e_z
        e_final[i] = ω/2pi * e_t
    end

    return posz,e_final
end

# path      = raw"D:\DualSPHysics_v5.0\examples\main\18_WavesTest\08_DpSPS\01_GenericWaveTank_mDBC_1.2_Dp-1\energyflux"
# pres_file = raw"EF10pts_01_GenericWaveTank_mDBC_1_Press.csv"
# vel_file  = raw"EF10pts_01_GenericWaveTank_mDBC_1_Vel.csv"
# fname_pres     = joinpath(path,pres_file)
# fname_vel     = joinpath(path,vel_file)

# P_posx,P_posz,P_header,P_df = read_pres(fname_pres)
# V_pos,V_header_x,V_df_x     = read_vel(fname_vel)
# posz,e_final = calc_energy(P_posx,P_posz,P_header,P_df,V_df_x,ids)
# plot(posz,e_final[1])
# ## Working on height variance

## Make small interface
function get_filename_gui()
    win = GtkWindow("Save File",300,50)

    ent = GtkEntry()
    set_gtk_property!(ent,:text,"Insert Filename")

    push!(win,ent)
    Gtk.showall(win)

    str = ""
    if isinteractive()
        c = Condition()
        signal_connect(win, "key-press-event") do widget, event
            if event.keyval == 65293 #Enter
                str = get_gtk_property(ent,:text,String)
                notify(c)
            end
        end
        wait(c)
    end

    Gtk.destroy(win)

    return str
end

parent_dir = raw"D:\DualSPHysics_v5.0\examples\main\19_FinalInitial"
dirs = filter(x -> isdir(joinpath(parent_dir, x)), readdir(parent_dir))

folders = joinpath.(joinpath.(parent_dir,dirs),"Simulation_Dp0.01")

#folders = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER,select_multiple=true)
#df = make_plots(folders)
filename = get_filename_gui()
plot_df_res(folders,filename)


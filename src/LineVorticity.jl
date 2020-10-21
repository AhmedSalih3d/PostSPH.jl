using Statistics #to use "mean"
using Plots
using Plots.PlotMeasures #for "margin"
using DataFrames
using CSV
using Gtk
using PGFPlotsX
using PlotThemes
using Colors
pgfplotsx()

theme(:vibrant,
             fglegend = colorant"#000000", #legend
             guidefontcolor = colorant"#000000",
             showaxis = true,
             #formatter = :scientific,
             aspect_ratio = :equal, #:equal
             #widen=true,
             #ylims=(0.1,100),xlims=(0,100),
             thickness_scaling = 1.25, 
             framestyle=:box,
             size=(1300,400),
             xtickfontsize=22,
             ytickfontsize=22,
             yguidefontsize=22,
             xguidefontsize=22,
             legendfontsize=22,
             titlefontsize=22,
             linewidth = 2.5,
             margin=10mm,
               )

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


searchdir(path,key) = filter(x->occursin(key,x), readdir(path))

function read_vor(top_path)

    sub_path = "vorticity"
    file_ext = "Vor.csv"

    f_path   = joinpath(top_path,sub_path)

    csv_file = searchdir(f_path,file_ext)[1]

    fname = joinpath(f_path,csv_file)
    file = open(fname,"r")

    f_parse_f64(f) = parse.(Float64,split(f,";")[3:end])

    f_pos    = readline(file)
    f_header = readline(file)

    f_pos_parsed = f_parse_f64(f_pos)

    vor_pos = reshape(f_pos_parsed,(3,Int(length(f_pos_parsed)/3)))

    vor_header = replace.(split(f_header,";"),Ref(r" \[m\]"=>""))
    vor_header = replace.(vor_header,Ref(r" \[s\]"=>""))
    vor_header = replace.(vor_header,Ref(r" \[1/s\]"=>""))

    vor_header_y = findall( x -> occursin("y", x), vor_header)

    vor_df = CSV.read(fname,DataFrame; delim=";",skipto=3,header=false)
    rename!(vor_df,Symbol.(vor_header))

    time   = vor_df[!,2]

    vor_df_y = vor_df[!,vor_header_y]

    close(file)
    return time,vor_pos,vor_header_y,vor_df_y
end

function do_plots(folders)
    all_plots = empty([], Plots.Plot)

    for folder in folders
        sub_path = "vorticity"

        time, vor_pos,vor_header,vor_df_y = read_vor(folder)

        #Abs due to magnitude and only one component
        p = plot(vor_pos[1,:],abs.(convert(Array,eachrow(vor_df_y)[1])),label="0",legendtitle="Time [s]",legend=:topright)
        title!(basename(folder))
        xlabel!("Horizontal Position [m]")
        isodd(length(all_plots))
        ylabel!("Vorticity [1/s]")
        xlims!(-13,7)
        ylims!(0,45)
        ids = [201 401 601 801]
        for i in ids
            plot!(vor_pos[1,:],abs.(convert(Array,eachrow(vor_df_y)[i])), label="$(round(time[i];digits=0))")
        end
        display(p)

        save_path = joinpath( joinpath(folder,sub_path) , "VorticityBreaking.pdf")
        savefig(save_path)

        push!(all_plots,p)
    end
    return all_plots
end

folders = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER,select_multiple=true)
all_plots = do_plots(folders)
n_plots   = length(all_plots)
if isodd(n_plots)
    n_plots = n_plots + 1
    push!(all_plots,plot(framestyle = :none))
end
ev_id = 2:2:n_plots
un_id = 1:2:n_plots
xlabel!.(all_plots[1:(n_plots - 2)],"")
ylabel!.(all_plots[ev_id],"")


sz = Int(n_plots/2)
l = @layout grid(sz,2)


pall = plot(all_plots..., layout = l,size=(800*sz,800*sz))
savefig(pall,"Test.pdf")

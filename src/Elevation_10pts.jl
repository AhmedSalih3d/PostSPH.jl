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


searchdir(path,key) = filter(x->occursin(key,x), readdir(path))


function plot_it(path,fname,df,a_left,a_right,e_left,e_right,posx,col_mean,mh,mht)
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
    group_size="1 by 5",
    #yticklabels_at="edge left",
    #xticklabels_at="edge bottom",
    vertical_sep="45pt",
    },
    no_markers
    },
    {legend_pos  = "outer north east", ylabel = "Amplitude [m]", height = "8.25cm", width = "15cm",grid,xticklabel=raw"\empty",ytick = [-0.3,-0.15,0,0.15,0.3],legend_style="{name=leg1,draw=none}"},
    hlines,
    PlotInc({color = cc[2]}, Table(df[!,"Time"], df[!,"Elevation_0"])),
    LegendEntry(raw"Elevation\_0"),
    PlotInc({color = cc[3]},Table(df[!,"Time"], df[!,"Elevation_1"])),
    LegendEntry(raw"Elevation\_1"),
    PlotInc({color = cc[4]},Table(df[!,"Time"], df[!,"Elevation_2"])),
    LegendEntry(raw"Elevation\_2"),
    PlotInc({color = cc[5]},Table(df[!,"Time"], df[!,"Elevation_3"])),
    LegendEntry(raw"Elevation\_3"),
    PlotInc({color = cc[6]},Table(df[!,"Time"], df[!,"Elevation_4"])),
    LegendEntry(raw"Elevation\_4"),
    PlotInc({color = cc[7],solid},Table(df[!,"Time"], df[!,"Elevation_5"])),
    LegendEntry(raw"Elevation\_5"),
    PlotInc({color = cc[1], thick, solid},Table(df[!,"Time"], a_left)),
    LegendEntry(raw"Mean"),
    {legend_pos  = "outer north east",xlabel = "Time [s]", ylabel = "Amplitude [m]", height = "8.25cm", width = "15cm",grid,ytick = [-0.3,-0.15,0,0.15,0.3],legend_style="{name=leg2,draw=none}"},
    hlines,
    PlotInc({color = cc[2]},Table(df[!,"Time"], df[!,"Elevation_6"])),
    LegendEntry(raw"Elevation\_6"),
    PlotInc({color = cc[3]},Table(df[!,"Time"], df[!,"Elevation_7"])),
    LegendEntry(raw"Elevation\_7"),
    PlotInc({color = cc[4]},Table(df[!,"Time"], df[!,"Elevation_8"])),
    LegendEntry(raw"Elevation\_8"),
    PlotInc({color = cc[5]},Table(df[!,"Time"], df[!,"Elevation_9"])),
    LegendEntry(raw"Elevation\_9"),
    PlotInc({color = cc[6]},Table(df[!,"Time"], df[!,"Elevation_10"])),
    LegendEntry(raw"Elevation\_10"),
    PlotInc({color = cc[7],solid},Table(df[!,"Time"], df[!,"Elevation_11"])),
    LegendEntry(raw"Elevation\_11"),
    PlotInc({color = cc[1], thick, solid},Table(df[!,"Time"], a_right)),
    LegendEntry(raw"Mean"),
    #Height Mean
    {legend_pos  = "outer north east",title=raw"Height Mean - $\langle h \rangle_t$", xlabel = "Elevation X", ylabel = raw"Height Mean - $\left[\SI{}{\meter}\right]$", height = "8.25cm", width = "15cm",grid,legend_style="{name=leg5,draw=none}",xtick=posx,yshift=-25},
    PlotInc({color = cc[1],"ultra thick"},Table(posx, col_mean)),
    LegendEntry(raw"Height Mean"),
    #Height variance
    {legend_pos  = "outer north east",title=raw"Height Variance - $\langle \left( h - \langle h \rangle \right)^2 \rangle_t$", xlabel = "Elevation X", ylabel = raw"Height Variance - $\left[\SI{}{\meter\squared}\right]$", height = "8.25cm", width = "15cm",grid,legend_style="{name=leg6,draw=none}",xtick=posx,yshift=-25},
    PlotInc({color = cc[1],"ultra thick"},Table(posx, mht)),
    LegendEntry(raw"Height Variance"),
    #Energy plot
    {legend_pos  = "outer north east", title=raw"Energy Mean Calculation -$\frac{1}{2}\rho g a^2$",xlabel = "Time [s]", ylabel = "Energy [J]", height = "7.5cm", width = "15cm",grid,legend_style="{name=leg7,draw=none}",yshift=-25},
    HLine({color = "gray", style ="{dashed}"}, 0),
    PlotInc({color = cc[1],thick},Table(df[!,"Time"], e_left)),
    LegendEntry(raw"Energy Left Edge"),
    PlotInc({color = cc[3],thick},Table(df[!,"Time"], e_right)),
    LegendEntry(raw"Energy Right Edge"),
    )

    push!(tp,p1)


    push!(tp,table_tex(path,posx))

    pgfsave(joinpath(path,fname)*".pdf", tp; include_preamble = true, dpi = 600)

    println("Done with $fname")
end


function make_plots(top_path)
    fpath = "elevation"

    final_path = joinpath.(top_path,fpath)

    println("Remember you have fixed rho and g!")
    rho = 1000
    g   = 9.8065

    file_ext = "csv"
    #file_ext = "grid"
    for path in final_path
        csv_file = searchdir(path,file_ext)

        println(csv_file[1])

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
        n_each = Int(n_col/2)

        left_range  = n_start:(n_ign + n_each)
        right_range = (left_range[end]+1):(left_range[end]+n_each)

        df = CSV.read(true_path,DataFrame; delim=";",skipto=5,header=false)
        rename!(df,Symbol.(header))

        mean_cols(df,range_val) = dropdims(mean(convert(Matrix,df[!,range_val]),dims=2),dims=2)

        a_left  = mean_cols(df,left_range)
        a_right = mean_cols(df,right_range)

        e_left  = 0.5 * rho * g * a_left.^2
        e_right = 0.5 * rho * g * a_right.^2

        col_mean = mean.(eachcol(df[:,n_start:end]))
        mh  = similar(col_mean)
        mht = similar(col_mean)
        for i = 1:length(col_mean)
            mh[i]  = mean((df[!,n_ign+i] .- col_mean[i]))
            mht[i] = mean((df[!,n_ign+i] .- col_mean[i]).^2)
        end

        plot_it(path,csv_file[1],df,a_left,a_right,e_left,e_right,posx,col_mean,mh,mht)

    end
end

folders = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER,select_multiple=true)
make_plots(folders)


## Working on height variance

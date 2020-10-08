using Statistics #to use "mean"
using Plots
using PlotThemes
using DataFrames
using CSV
using Gtk
gr()



energy_path = "energy_plots"

file_path = "df.csv"

sma(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]

function plot_test(final_path,sym)
    p = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
        for path in final_path
            data = CSV.read(path, DataFrame;)
            plot!(data.t,data[!,sym],label=splitpath(dirname(path))[end-1])
        end
    display(p)
end

function plot_test_tot_1(final_path)
    p = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
        for path in final_path
            data = CSV.read(path, DataFrame;)
            plot!(data.t,data.ke1 .+ data.pe1,label=splitpath(dirname(path))[end-1])
        end
    display(p)
end



function plot_test_sma(final_path,sym)
    p = plot(legend=:outerbottom,foreground_color_legend = nothing,background_color_legend = nothing)
        for path in final_path
            data = CSV.read(path, DataFrame;)
            plot!(sma(data[!,sym],100),label=splitpath(dirname(path))[end-1])
        end
    display(p)
end

top_path = open_dialog("Select Dataset Folder", action=GtkFileChooserAction.SELECT_FOLDER,select_multiple=true)
final_path = joinpath.(joinpath.(top_path,energy_path),file_path)
title!("Right Edge KE - ARTIFICIAL")
ylabel!("Joule")
plot_test_tot_1(final_path)
#plot_test_sma(final_path,:de2)

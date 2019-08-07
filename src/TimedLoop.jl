using Statistics
using Plots

cd("C:\\Users\\ahmed\\Documents\\DualSPHysics4.4\\DualSPHysics_v4.4\\examples\\main\\03_MovingSquare\\CaseMovingSquare_out")
Bodies = MkArray("CaseMovingSquare.xml")
MovingSquare = Bodies[2][1]
cd("C:\\Users\\ahmed\\Documents\\DualSPHysics4.4\\DualSPHysics_v4.4\\examples\\main\\03_MovingSquare\\CaseMovingSquare_out\\data")

function test()
    plot(xlabel="Time [s]",ylabel="x-position [m]")
    i = 1
    n = 1
    while true
        if i == 1
            difFiles = _dirFiles()
            i = 0
        else
            curFiles = _dirFiles()
            sleep(2.5)
            newFiles = _dirFiles()
            difFiles = setdiff(newFiles,curFiles)
        end
        data     = readBi4Body(MovingSquare,Points,difFiles)
        time     = readBi4Time(difFiles)
        y        = getindex.(mean.(data),1)[:]
        pp = plot!(time,y,label="")
        display(pp)
    end
end

#ReadBi4Rhop.jl
cd("C:/Users/ahmed/Documents/DualSPHysics4.4/DualSPHysics_v4.4/examples/main/03_MovingSquare/CaseMovingSquare_out/data")


function bi4test()
    ft = open("Part_0000.bi4",read=true)

    readuntil(ft,"Rhop")
    readuntil(ft,"Rhop")

    for i = 1:2
        read(ft,Int32)
    end

    n = read(ft,Int32)

    read(ft,Int32)

    k = zeros(Float32,n)
    for i = 1:n
        global k[i] = read(ft,Float32)
    end
    return k
end

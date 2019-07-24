#if typ == Idp || typ == Rhop
#    j[i] = zeros(catTypeBi4[typ], n)
#    _transferDataBi4(ft,j[i])
#elseif typ == Points || typ == Vel
#    j[i] = zeros(catTypeBi4[typ], (n,catColBi4[typ]))
#    _transferDataBi4(ft,j[i])
#end

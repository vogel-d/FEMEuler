function l2g(degF::degF{N} where N,fid::Int64)
    #globale Nummerierung der vertices einer Entität aus D besteht aus
    #beliebigen Zahlen, startend "unten links", gegen den Uhrzeigersinn velaufend
    #vertices der Entität fid aus D werden festgestellt und zurückgegeben,
    #durch (notwendige) korrekte Eingabe der Vertices innerhalb einer Entität aus D
    #ist dies bereits sortiert wie die lokale Nummerierung
    #(also gegen den Uhrzeigersinn, startend "unten links")
    return degF.incidence[degF.offset[fid]:degF.offset[fid+1]-1];
end

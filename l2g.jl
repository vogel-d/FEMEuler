function l2g(degF::degF{N,S} where N where S,c::Int)
    #globale Nummerierung der vertices einer Entität aus D besteht aus
    #beliebigen Zahlen, startend "unten links", gegen den Uhrzeigersinn velaufend
    #vertices der Entität c aus D werden festgestellt und zurückgegeben,
    #durch (notwendige) korrekte Eingabe der Vertices innerhalb einer Entität aus D
    #ist dies bereits sortiert wie die lokale Nummerierung
    #(also gegen den Uhrzeigersinn, startend "unten links")
    return degF.incidence[degF.offset[c]:degF.offset[c+1]-1];
end


function l2g!(globalNum::Array{Int,1},degF::degF{N,S} where N where S,c::Int)
    #globale Nummerierung der vertices einer Entität aus D besteht aus
    #beliebigen Zahlen, startend "unten links", gegen den Uhrzeigersinn velaufend
    #vertices der Entität c aus D werden festgestellt und zurückgegeben,
    #durch (notwendige) korrekte Eingabe der Vertices innerhalb einer Entität aus D
    #ist dies bereits sortiert wie die lokale Nummerierung
    #(also gegen den Uhrzeigersinn, startend "unten links")
    z=1;
    for i=degF.offset[c]:(degF.offset[c+1]-1)
        globalNum[z]=degF.incidence[i];
        z+=1;
    end
    return nothing;
end

function l2g!(globalNum::Array{Int,1},degF::degF{N,S} where N where S,cells::Array{Int,1})
    #globale Nummerierung der vertices einer Entität aus D besteht aus
    #beliebigen Zahlen, startend "unten links", gegen den Uhrzeigersinn velaufend
    #vertices der Entität c aus D werden festgestellt und zurückgegeben,
    #durch (notwendige) korrekte Eingabe der Vertices innerhalb einer Entität aus D
    #ist dies bereits sortiert wie die lokale Nummerierung
    #(also gegen den Uhrzeigersinn, startend "unten links")
    z=1;
    for c in cells
        for i in degF.offset[c]:(degF.offset[c+1]-1)
            globalNum[z]=degF.incidence[i];
            z+=1;
        end
    end
    return nothing;
end

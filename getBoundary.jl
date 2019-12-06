function getBoundary(m::mesh)
  #boundary = alle Entitäten aus D-1, die nur zu einer Entität aus D inzident sind
    dim=m.topology.D;
    #"D-1D" wird erzeugt, um Nachbarn festzustellen
    meshConnectivity!(m,dim-1,dim);
    off=m.topology.offset["$(dim-1)$dim"];
    inc=m.topology.incidence["$(dim-1)$dim"];
    b=Set{Int64}();
    for e in 1:(length(off)-1)
        if off[e+1]-off[e]==1
            #immer, wenn die Grenzen im Offset sich nur um 1 unterscheiden
            #wird die Entitäts-ID in boundary gespeichert
            push!(b,e)
        end
    end

    return b
end

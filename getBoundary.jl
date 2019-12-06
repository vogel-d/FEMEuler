function getBoundary(m::mesh)
  #boundary = alle Kanten, die nur zu einer Fläche inzident sind
    meshConnectivity!(m,1,2);
    off=m.topology.offset["12"];
    inc=m.topology.incidence["12"];
    offv=m.topology.offset["10"];
    incv=m.topology.incidence["10"];
    be=Set{Int}();
    bv=Set{Int}();
    for e in 1:m.topology.size[2]
        if off[e+1]-off[e]==1
            #immer, wenn die Grenzen im Offset sich nur um 1 unterscheiden
            #wird die Entitäts-ID in boundary gespeichert
            v=incv[offv[e]:offv[e+1]-1];
            push!(be,e)
            union!(bv,v)
        end
    end

    return be, bv
end

function getBoundary(m::mesh)
  #boundary = alle Entitäten aus D-1, die nur zu einer Entität aus D inzident sind
    dim=m.topology.D;
    #"D-1D" wird erzeugt, um Nachbarn festzustellen
    meshConnectivity!(m,dim-1,dim);
    off=m.topology.offset["$(dim-1)$dim"];
    inc=m.topology.incidence["$(dim-1)$dim"];

    b=Dict{Int64, Array{Int64,1}}();
    for i in 1:(length(off)-1)
        if off[i+1]-off[i]==1
            #immer, wenn die Grenzen im Offset sich nur um 1 unterscheiden
            #wird die Entitäts-ID in boundary gespeichert
            fi=inc[off[i]];
            if haskey(b, fi)
                push!(b[fi],i);
            else
                b[fi]=[i];
            end
        end
    end

    return b
end

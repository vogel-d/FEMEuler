function getStencil(m::mesh, order::Int)
    stencil=Array{Array{Int,1},1}(undef,m.topology.size[3]);
    meshConnectivity!(m,2,1)
    meshConnectivity!(m,1,2)
    ince=m.topology.incidence["21"]
    offe=m.topology.offset["21"]
    incf=m.topology.incidence["12"]
    offf=m.topology.offset["12"]
    edges=Set{Int}();
    ch=Int[]
    cells=Set{Int}();
    for f in 1:m.topology.size[3]
        cells=Set{Int}(f);
        ocells=Int[f];
        for i in 1:order
            for c in ocells
                for e in ince[offe[c]:offe[c+1]-1]
                    if !in(e,edges)
                        ce=incf[offf[e]:offf[e+1]-1]
                        if length(ce)==2
                            if c==ce[1]
                                push!(cells,ce[2])
                                push!(ch,ce[2])
                            elseif c==ce[2]
                                push!(cells,ce[1])
                                push!(ch,ce[1])
                            end
                        end
                        push!(edges,e)
                    end
                end
            end
            ocells=ch;
            ch=Int[];
        end
        stencil[f]=collect(cells);
        edges=Int[];
    end
    return stencil
end

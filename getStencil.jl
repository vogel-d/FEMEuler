function getStencil(m::mesh, order::Int)
    stencil=Array{Array{Int,1},1}(undef,m.topology.size[3]);
    stencilBoundary=spzeros(Int,m.topology.size[3])
    meshConnectivity!(m,2,1)
    meshConnectivity!(m,1,2)
    ince=m.topology.incidence["21"]
    offe=m.topology.offset["21"]
    incf=m.topology.incidence["12"]
    offf=m.topology.offset["12"]
    incv=m.topology.incidence["10"]
    offv=m.topology.offset["10"]
    mcoord=m.geometry.coordinates;
    bE=copy(m.boundaryEdges);

    edges=Set{Int}();
    ch=Int[]
    cells=Set{Int}();
    boundcells=Set{Int}()
    boundind=Int[]
    ocells=Int[];
    for f in 1:m.topology.size[3]
        empty!(cells)
        push!(cells,f)
        empty!(ocells)
        push!(ocells,f)
        for i in 1:order
            for c in ocells
                for e in ince[offe[c]:offe[c+1]-1]
                    if !in(e,edges)
                        ce= @views incf[offf[e]:offf[e+1]-1]
                        if length(ce)==2
                            if c==ce[1]
                                push!(cells,ce[2])
                                push!(ch,ce[2])
                            elseif c==ce[2]
                                push!(cells,ce[1])
                                push!(ch,ce[1])
                            end
                        elseif length(ce)==1
                            if bE[e]<0
                                e2=-bE[e]
                                bE[e2]=-e
                                push!(edges,e2)
                                push!(cells,incf[offf[e2]])
                                push!(ch,incf[offf[e2]])
                                push!(boundcells,incf[offf[e2]])
                                coord=@views mcoord[:,incv[offv[e]:offv[e+1]-1]]
                                coord2=@views mcoord[:,incv[offv[e2]:offv[e2+1]-1]]
                                if coord[2,1]==coord[2,2]
                                    if coord[2,1]<coord2[2,1]
                                        stencilBoundary[c]=-1
                                    elseif coord[2,1]>coord2[2,1]
                                        stencilBoundary[c]=1
                                    end
                                elseif coord[1,1]==coord[1,2]
                                    if coord[1,1]<coord2[1,1]
                                        stencilBoundary[c]=-2
                                    elseif coord[1,1]>coord2[1,1]
                                        stencilBoundary[c]=2
                                    end
                                end
                            end
                        end
                        push!(edges,e)
                    end
                end
            end
            copy!(ocells,ch);
            empty!(ch);
        end
        stencil[f]=collect(cells);
        for c in boundcells
            boundind=findall(c,stencil[f])
            @. stencil[f][boundind] *= -1
        end
        empty!(boundcells)
        empty!(edges);
    end
    return stencil, stencilBoundary
end

function getStencil(m::mesh, order::Float64)
    stencil=Array{Array{Int,1},1}(undef,m.topology.size[3]);
    stencilBoundary=spzeros(Int,m.topology.size[3])
    meshConnectivity!(m,2,2)
    inc=m.topology.incidence["22"]
    off=m.topology.offset["22"]
    for f in 1:m.topology.size[3]
        stencil[f]=Int[]
        push!(stencil[f],f)
        append!(stencil[f],inc[off[f]:off[f+1]-1])
    end
    return stencil, stencilBoundary
end

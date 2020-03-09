function setOrientation!(m::mesh)
    meshConnectivity!(m,2,1)
    meshConnectivity!(m,1,2)

    ne=m.topology.size[m.topology.dim]
    nf=m.topology.size[m.topology.dim+1]

    visited=falses(ne);
    for e in 1:ne
        setOrientation!(m,visited,e,m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1])
    end

    @views ince=m.topology.incidence["21"]
    @views offe=m.topology.offset["21"]
    @views inc=m.topology.incidence["10"]
    @views off=m.topology.offset["10"]

    incfv=Int[]; incfe=Int[];
    for c in 1:nf
        edges=ince[offe[c]:offe[c+1]-1];
        v=0; vert=[0,0,0]; pedges=[0,0,0,0]; ie=1;
        while v==0 && ie<length(edges)
            v1=inc[off[edges[ie]]];
            v2=inc[off[edges[ie]+1]-1];
            for je in ie+1:length(edges)
                if v1==inc[off[edges[je]+1]-1]
                    v=v1;
                    vert[1]=inc[off[edges[ie]+1]-1]
                    vert[3]=inc[off[edges[je]]];
                    pedges[1]=edges[ie]
                    pedges[4]=edges[je]
                    re=setdiff(edges,pedges)
                    ve=inc[off[re[1]]:off[re[1]+1]-1]
                    se=(vert[1] .== ve)
                    if sum(se)==1
                        pedges[2]=re[1]
                        pedges[3]=re[2]
                        vert[2]=ve[findall(.!(se))[1]];
                    else
                        pedges[2]=re[2]
                        pedges[3]=re[1]
                        se=(vert[3] .!= ve)
                        vert[2]=ve[findall(se)[1]];
                    end
                    break;
                #=
                if v1==inc[off[edges[je]]]
                    v=v1;
                    vert[1]=inc[off[edges[ie]+1]-1]
                    vert[3]=inc[off[edges[je]+1]-1];
                    pedges[1]=edges[ie]
                    pedges[4]=edges[je]
                    re=setdiff(edges,pedges)
                    ve=inc[off[re[1]]:off[re[1]+1]-1]
                    se=(vert[1] .== ve)
                    if sum(se)==1
                        pedges[2]=re[1]
                        pedges[3]=re[2]
                        vert[2]=ve[findall(.!(se))[1]];
                    else
                        pedges[2]=re[2]
                        pedges[3]=re[1]
                        se=(vert[3] .!= ve)
                        vert[2]=ve[findall(se)[1]];
                    end
                    break;
                if v2==inc[off[edges[je]]]
                    v=v2;
                    vert[1]=inc[off[edges[je]+1]-1];
                    vert[3]=inc[off[edges[ie]]];
                    pedges[1]=edges[je]
                    pedges[4]=edges[ie]
                    re=setdiff(edges,pedges)
                    ve=inc[off[re[1]]:off[re[1]+1]-1]
                    se=(vert[1] .== ve)
                    if sum(se)==1
                        pedges[2]=re[1]
                        pedges[3]=re[2]
                        vert[2]=ve[findall(.!(se))[1]];
                    else
                        pedges[2]=re[2]
                        pedges[3]=re[1]
                        se=(vert[3] .!= ve)
                        vert[2]=ve[findall(se)[1]];
                    end
                    break;
                =#
                else
                    continue;
                end
            end
            ie+=1;
        end
        append!(incfe,pedges);
        pushfirst!(vert,v);
        coord=m.geometry.coordinates[:,vert]
        e0=coord[:,2]-coord[:,1];
        e1=coord[:,4]-coord[:,1];
        n1=cross(e0,e1)
        s=sign(dot(n1,coord[:,1]))
        #n=transformation(m,coord,0.5,0.5)
        append!(incfv,vert);
        push!(m.orientation,s)
    end
    m.topology.incidence["20"]=incfv;
    m.topology.incidence["21"]=incfe;
    return nothing
end

function setOrientation!(m::mesh,visited::BitArray{1},e::Int,v::Array{Int,1})
    if visited[e]
        if m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1]!=v
            error("Möbius strip found.")
        end
    else
        @views cells=m.topology.incidence["12"][m.topology.offset["12"][e]:m.topology.offset["12"][e+1]-1]
        @views ince=m.topology.incidence["21"]
        @views offe=m.topology.offset["21"]
        @views incv=m.topology.incidence["20"]
        @views offv=m.topology.offset["20"]

        visited[e]=true
        m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1]=v;
        for c in cells
            edges=ince[offe[c]:offe[c+1]-1];
            vert=incv[offv[c]:offv[c+1]-1][[1,2,4,3]];

            #Evtl. Auslagern in Funktion für allg. Behandlung?
            ie=findall(e,edges)[1];
            isodd(ie) ? e2=edges[ie+1] : e2=edges[ie-1];

            posv=findall(v[1],vert)
            append!(posv,findall(v[2],vert))

            v2=setdiff([1,2,3,4],posv)
            if posv[1]<posv[2]
                v2=vert[v2];
            else
                v2=vert[v2[[2,1]]];
            end

            !visited[e2] && setOrientation!(m,visited,e2,v2)
        end
    end
    return nothing
end

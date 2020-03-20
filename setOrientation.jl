function setOrientation!(m::mesh)
    meshConnectivity!(m,2,1)
    meshConnectivity!(m,1,2)

    ne=m.topology.size[m.topology.dim]
    nf=m.topology.size[m.topology.dim+1]

    incef=zeros(Int,2*ne);
    visited=falses(ne);
    for e in 1:ne
        setOrientation!(m,visited,incef,e,m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1])
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
    m.topology.incidence["12"]=incef;
    return nothing
end

function setOrientation!(m::mesh,visited::BitArray{1},inc::Array{Int,1},e::Int,v::Array{Int,1})
    if visited[e]
        if m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1]!=v
            error("Möbius strip found.")
        end
    else
        @views offc=m.topology.offset["12"]
        @views incc=m.topology.incidence["12"]
        @views ince=m.topology.incidence["21"]
        @views offe=m.topology.offset["21"]
        @views incv=m.topology.incidence["20"]
        @views offv=m.topology.offset["20"]

        visited[e]=true
        m.topology.incidence["10"][m.topology.offset["10"][e]:m.topology.offset["10"][e+1]-1]=v;
        permc=(2,1);
        if iszero(inc[offc[e]])
            cells=incc[offc[e]:offc[e+1]-1];
            inc[offc[e]:offc[e+1]-1]=cells;
        else
            cells=inc[offc[e]:offc[e+1]-1];
        end
        for ic in 1:length(cells)
            c=cells[ic];
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

            cells2=incc[offc[e2]:offc[e2+1]-1];
            if c==cells2[1]
                inc[offc[e2]:offc[e2+1]-1]=cells2[[permc[ic],ic]]
            elseif c==cells2[2]
                inc[offc[e2]:offc[e2+1]-1]=cells2[[ic,permc[ic]]]
            end
            !visited[e2] && setOrientation!(m,visited,inc,e2,v2)
        end
    end
    return nothing
end

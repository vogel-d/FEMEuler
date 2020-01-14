
function plotFEM(m::mesh, key::Symbol, showann::Bool=true; showvertices::Bool=showann,
    sizevertices::Float64=2.5, linecolor::Symbol=:blue,
    colorann::Array{Symbol,1}=[:grey, :darkgrey, :black],
    positionann::Array{Symbol,1}=[:bottom, :bottom, :auto],
    sizeann::Array{Int64,1}=[8,8,8])

    ordEdgesB, nebP, nebC=getOrderBoundary(m.boundaryEdges);
    ordVerticesB, nvbP, nvbC=getOrderBoundary(m.boundaryVertices);
    if m.meshType==4
        phi, divphi, gradphi, cm, refFace, refEdge, refVert=getQuadElementProperties(key);
    elseif m.meshType==3
        phi, divphi, gradphi, cm, refFace, refEdge, refVert=getTriElementProperties(key);
    end

    ince=m.topology.incidence["10"];

    incf=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    coord=m.geometry.coordinates;
    nv, ne, nf=m.topology.size[1:3];
    x=Array{Float64,2}(undef,2,ne);
    y=Array{Float64,2}(undef,2,ne);
    coorde=Array{Float64,2}(undef,2,ne);

    for k in 1:ne
        i=ince[[(2*k-1),2*k]];
        x[:,k]=coord[1,i];
        y[:,k]=coord[2,i];
        coorde[1,k]=0.5*sum(x[:,k]);
        coorde[2,k]=0.5*sum(y[:,k]);
    end

    coordf=Array{Float64,2}(undef,2,nf);
    ng=off[2]-off[1];

    for k in 1:nf
        i=incf[off[k]:(off[k+1]-1)];
        coordf[1,k]=sum(coord[1,i])/ng;
        coordf[2,k]=sum(coord[2,i])/ng;
    end

    p=plot(x,y, c=linecolor, linewidth=0.5, legend=false, xlabel="x in m", ylabel="z in m");
    if showann
        t=[];
        o=[1];
        te=[];
        oe=[1];
        tf=[];
        of=[1];
        for v in 1:nv
            if m.boundaryVertices[v]>=0
                for d in 1:refVert
                    push!(t,nf*refFace+(ne-nebP)*refEdge+refVert*(ordVerticesB[v]-1)+d)
                end
            elseif m.boundaryVertices[v]<0
                for d in 1:refVert
                    push!(t,nf*refFace+(ne-nebP)*refEdge+refVert*(ordVerticesB[-m.boundaryVertices[v]]-1)+d)
                end
            end
            push!(o,length(t)+1);
        end
        for e in 1:ne
            if m.boundaryEdges[e]>=0
                for d in 1:refEdge
                    push!(te,nf*refFace+refEdge*(ordEdgesB[e]-1)+d)
                end
            elseif m.boundaryEdges[e]<0
                for d in 1:refEdge
                    push!(te,nf*refFace+refEdge*(ordEdgesB[-m.boundaryEdges[e]]-1)+d)
                end
            end
            push!(oe,length(te)+1);
        end
        for f in 1:nf
            for d in 1:refFace
                push!(tf, refFace*(f-1)+d)
            end
            push!(of,length(tf)+1);
        end
        if !isempty(t)
            a=[(coord[1,i], coord[2,i], text("$(t[o[i]:o[i+1]-1])", colorann[1], positionann[1], sizeann[1])) for i in 1:nv];
        else
            a=[(coord[1,i], coord[2,i], text("$i", colorann[1], positionann[1], sizeann[1])) for i in 1:nv];
        end
        if !isempty(te)
            ae=[(coorde[1,i], coorde[2,i], text("$(te[oe[i]:oe[i+1]-1])", colorann[2], positionann[2], sizeann[2])) for i in 1:ne];
            append!(a,ae)
        end
        if !isempty(tf)
            af=[(coordf[1,i], coordf[2,i], text("$(tf[of[i]:of[i+1]-1])", colorann[3], positionann[3], sizeann[3])) for i in 1:nf];
            append!(a,af);
        end
        annotate!(a);
    end

    if showvertices
        plot!(coord[1,:], coord[2,:], seriestype=:scatter, c=:black, markersize=sizevertices);
    end

    return p
end

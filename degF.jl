struct degF{N,Space}
    # size: in assembLoad, assembMassRho, advectionStiff, assembStiff, applyStartValues, embed, testFunctions, projectPressure, projectRecovery
    # in generateBoundary.jl, plotFEM.jl
    numB::Int;
    num::Int;

    incidence::Array{Int,1}; # in l2g.jl, correctVelocity.jl, plotSolution.jl, plotSolutionGif.jl #notwendig
    offset::Array{Int,1}; # in l2g.jl, correctVelocity.jl, plotSolution.jl, plotSolutionGif.jl #evtl. weglassen und Schrittweite speichern? ->Dann nicht für Gitter mit gemischten Elementen, z.B. Dreiecke und Rechtecke

    phi::Array{Array{Float64,2},N};
    divphi::Array{Array{Float64,2},1};
    gradphi::Array{Array{Float64,2},2};
end


function degF(m::mesh, femType::Symbol, ordEdgesB::Array{Int,1}, nebP::Int, nebC::Int,
        ordVerticesB::Array{Int,1}, nvbP::Int, nvbC::Int, kubPoints::Array{Float64,2})
    nf=m.topology.size[3];
    ne=m.topology.size[2];
    nv=m.topology.size[1];
    meshConnectivity!(m,2,1);
    incfe=m.topology.incidence["21"];
    offfe=m.topology.offset["21"];
    incfv=m.topology.incidence["20"];
    offfv=m.topology.offset["20"];
    incev=m.topology.incidence["10"];
    offev=m.topology.offset["10"];
    #Gerade noch nur für gleichmäßige Meshes
    #Verallgemeinern durch Rauskürzen von nef und nvf und variablen Erstellen von off & Erweitern von getElementProperties
    nef=offfe[2]-offfe[1];
    nvf=offfv[2]-offfv[1];
    phi, divphi,  gradphi, refFace, refEdge, refVert=getElementProperties(femType, kubPoints, m.meshType, m.geometry.dim)

    ndegF=refFace+nef*refEdge+nvf*refVert
    inc=zeros(Int, nf*ndegF);
    off=collect(1:ndegF:nf*ndegF+1);
    for f in 1:nf
        for d in 1:refFace
            inc[off[f]+d-1]=refFace*(f-1)+d;
        end

        vert=incfv[offfv[f]:offfv[f+1]-1];
        zv=off[f]+refFace+nef*refEdge
        for v in vert
            if m.boundaryVertices[v]>=0
                for d in 1:refVert
                    inc[zv]=nf*refFace+(ne-nebP)*refEdge+refVert*(ordVerticesB[v]-1)+d; #evtl in zwei Konstanten speichern
                    zv+=1;
                end
            elseif m.boundaryVertices[v]<0
                for d in 1:refVert
                    inc[zv]=nf*refFace+(ne-nebP)*refEdge+refVert*(ordVerticesB[-m.boundaryVertices[v]]-1)+d; #evtl in zwei Konstanten speichern
                    zv+=1;
                end
            end
        end

        edges=incfe[offfe[f]:offfe[f+1]-1];

        if m.topology.dim==m.geometry.dim
            #Sortieren von edges in  richtige Reigenforge nach Knoten
            push!(vert,vert[1]);
            ind=zeros(Int,length(edges)); #ALLOCATION
            for i in 1:length(edges)
                ve=incev[offev[edges[i]]:offev[edges[i]+1]-1];
                for j in 1:length(vert)-1
                    if ve[1]==vert[j]
                        if ve[2]==vert[j+1]
                            ind[j]=i;
                        else
                            continue;
                        end
                    elseif ve[1]==vert[j+1]
                        if ve[2]==vert[j]
                            ind[j]=i;
                        else
                            continue;
                        end
                    end
                end
            end
            edges=edges[ind];
        end

        ze=off[f]+refFace;
        for e in edges
            if m.boundaryEdges[e]>=0
                for d in 1:refEdge
                    inc[ze]=nf*refFace+refEdge*(ordEdgesB[e]-1)+d; #evtl in zwei Konstanten speichern
                    ze+=1;
                end
            elseif m.boundaryEdges[e]<0
                for d in 1:refEdge
                    inc[ze]=nf*refFace+refEdge*(ordEdgesB[-m.boundaryEdges[e]]-1)+d; #evtl in zwei Konstanten speichern
                    ze+=1;
                end
            end
        end
    end
    nb=nf*refFace+(ne-nebP)*refEdge+(nv-nvbP)*refVert;
    n=nb-nebC*refEdge-nvbC*refVert;
    degF{ndims(phi),getSpace(femType)}(nb,n, inc, off, phi, divphi, gradphi);
end

#TODO: DELETE GETSPACE.jl

function getSpace(femType::Symbol)

    H1=Set([:DG0,:DG1,:DG2,:P1,:P2]);
    H1div=Set([:RT0,:RT1,:RT0B,:RT1B]);
    H1xH1=Set([:VecDG1,:VecP1,:VecDG1S,:VecP1S]);

    if in(femType,H1)
        return :H1
    elseif in(femType,H1div)
        return :H1div
    elseif in(femType,H1xH1)
        return :H1xH1
    else
        error("Bitte für $femType Raum spezifizieren.")
    end
end

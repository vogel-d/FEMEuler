mutable struct degF{N}
    # size: in assembLoad, assembMassRho, advectionStiff, assembStiff, applyStartValues, embed, testFunctions, projectPressure, projectRecovery
    # in generateBoundary.jl, plotFEM.jl
    numB::Int;
    num::Int;

    #edgeBoundaryIncidence::SparseVector{Array{Int,1},Int}; #in generateBoundary.jl
    #boundaryEdgeIncidence::Array{Array{Int,1},1}; #in generateBoundary.jl

    incidence::Array{Int,1}; # in l2g.jl, correctVelocity.jl, plotSolution.jl, plotSolutionGif.jl #notwendig
    offset::Array{Int,1}; # in l2g.jl, correctVelocity.jl, plotSolution.jl, plotSolutionGif.jl #evtl. weglassen und Schrittweite speichern?

    #referenceBoundary::Array{Float64,2}; #in setEdgeData.jl ->setEdgeData in degF/femProblem, dann Speichern nicht notwendig oder extra Funktion, die es abhängig von FE-Raum aussucht
    #-> cm in getElementProperties
    phi::Array{Array{Float,2},N};
    divphi::Array{Array{Float,2},1};
    gradphi::Array{Array{Float,2},2};
    components::Array{Int,1}; #in plotSolution.jl und vtk.jl
end


function degF(m::mesh, femType::Symbol, boundaryEdges::Set{Int}, boundaryVertices::Set{Int}, kubPoints::Array{Float64,2})
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

    phi, divphi,  gradphi, comp,  refFace, refEdge, refVert=getQuadElementProperties(femType, kubPoints)
    ndegF=refFace+nef*refEdge+nvf*refVert
    inc=zeros(Int, nf*ndegF);
    off=collect(1:ndegF:nf*ndegF+1);

    #nbe=length(boundaryEdges); #Anzahl äußerer Kanten
    #nie=ne-nbe; #Anzahl innerer Kanten
    indEint=Array{Int,1}();
    indEbound=Array{Int,1}();
    indVint=Array{Int,1}();
    indVbound=Array{Int,1}();
    for f in 1:nf
        for d in 1:refFace
            inc[off[f]+d-1]=refFace*(f-1)+d;
        end

        vert=incfv[offfv[f]:offfv[f+1]-1];
        zv=off[f]+refFace+nef*refEdge
        for v in vert
            if !in(v,boundaryVertices)
                for d in 1:refVert
                    inc[zv]=nf*refFace+ne*refEdge+refVert*(v-1)+d; #evtl in zwei Konstanten speichern
                    push!(indVint,zv);
                    zv+=1;
                end
            else
                for d in 1:refVert
                    inc[zv]=nf*refFace+ne*refEdge+refVert*(v-1)+d; #evtl in zwei Konstanten speichern
                    push!(indVbound,zv);
                    zv+=1;
                end
            end
        end
        #Randbedingungen, bei Erstellen BoundaryEdges alle anliegenden Vert in anderes Set

        edges=incfe[offfe[f]:offfe[f+1]-1];
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
        ze=off[f]+refFace;
        for e in edges
            if !in(e,boundaryEdges)
                for d in 1:refEdge
                    inc[ze]=nf*refFace+refEdge*(e-1)+d; #evtl in zwei Konstanten speichern
                    push!(indEint,ze)
                    ze+=1;
                end
            else
                #Free Slip RB
                for d in 1:refEdge
                    inc[ze]=nf*refFace+refEdge*(e-1)+d; #evtl in zwei Konstanten speichern
                    push!(indEbound,ze)
                    ze+=1;
                end
            end
        end
    end

    #Sortieren der Randfreiheitsgrade nach hinten
    for i in 1:length(indEint)
        inc[indEint[i]]=
    end
    for i in 1:length(indEbound)
        inc[indEbound[i]]=
    end
    for i in 1:length(indVint)
        inc[indVint[i]]=
    end
    for i in 1:length(indVbound)
        inc[indVbound[i]]=
    end

    nb=nf*refFace+ne*refEdge+nv*refVert;
    n=nb-length(boundaryEdges)*refEdge;
    degF(nb,n, inc, off, phi, divphi, gradphi, comp);
end

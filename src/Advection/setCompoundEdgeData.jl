function setCompoundEdgeData!(p::femProblem, compVf::Symbol)
    m=p.mesh;
    degFVf=p.degFBoundary[p.femType[compVf][1]];
    mt=m.meshType;
    nVertices_CompoundElement=diff(m.topology.offset["20"][1:2])[1];
    refBound, edgeTypes=getCompoundElementProperties(p.femType[compVf][1],p.data.compoundData);
    #if mt==4
    #    edgeTypes=Dict([1,2]=>1,[2,3]=>2,[3,4]=>3,[1,4]=>4)
    #    coordref=[0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0]
    #else
    #    edgeTypes=Dict([1,2]=>1,[2,3]=>2,[1,3]=>3)
    #    coordref=[0.0 1.0 0.0; 0.0 0.0 1.0]
    #end
    clockwiseRotation = [0 1; -1 0];

    meshConnectivity!(m,1,2)
    offe=m.topology.offset["12"];
    ince=m.topology.incidence["12"];
    offv=m.topology.offset["10"];
    incv=m.topology.incidence["10"];
    offf=m.topology.offset["20"];
    incf=m.topology.incidence["20"];
    edgeNum=Int64[];
    edgeType=Int64[];
    cells=Int64[];
    globv=Int64[];
    off=Int64[1];
    periodic=Bool[];
    zo=1;
    for e in 1:m.topology.size[m.topology.dim]
        off1=offe[e];
        off2=offe[e+1];
        h1=false;
        switched=false;
        if (off2-off1)==1
            e1=m.boundaryEdges[e];
            if e1<0
                e1=-e1
                inc=[ince[offe[e1]], ince[off1]];
                h1=true;
            else
                continue;
            end
        else
            inc=(ince[off1:(off2-1)]);
        end

        #make sure global orientation is right ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #(so that global edge-normal points from inc[1] to inc[2])
        #(global edge-normal is the vector
        # from lower-number-vertex to higher-number-vertex
        # rotated 90° clockwise)
        # -> switch inc1 and inc2 if needed
        #TODO: think about how it works in 3d

        #get numbers of vertices
        edgeVertices = incv[offv[e]:offv[e]+1];
        #create vector pointing from vertex lower number to vertex higher number
        # (with respect to global edge orientation)
        sort!(edgeVertices);
        globalNormal = m.geometry.coordinates[:,edgeVertices[2]] .- m.geometry.coordinates[:,edgeVertices[1]]
        #rotate it clockwise (becomes global normal)
        globalNormal = clockwiseRotation*globalNormal;

        #compute centers of cells
        verticesInc1 = m.geometry.coordinates[:,incf[offf[inc[1]]:(offf[inc[1]+1]-1)]]
        verticesInc2 = m.geometry.coordinates[:,incf[offf[inc[2]]:(offf[inc[2]+1]-1)]]
        centerInc1 = 1/(nVertices_CompoundElement) * sum(verticesInc1,dims=2);
        centerInc2 = 1/(nVertices_CompoundElement) * sum(verticesInc2,dims=2);

        #define current normal pointing vom cell 1 to cell 2
        #old version
        #currentNormal = centerInc2 .- centerInc1;
        #note: if edge e is on periodic boundary (h1=true) the computed "currentNormal" will have
        #the opposite direction of the actual current normal (which is pointing out of the grid), so:
        #if h1
        #    currentNormal= -currentNormal;
        #end

        #new version (better handling at periodic boundary)
        #current normal is equivalent to a current normal pointing from the middle of the edge to cell 2
        centerEdge = 0.5* (m.geometry.coordinates[:,edgeVertices[2]] .+ m.geometry.coordinates[:,edgeVertices[1]]);
        currentNormal = centerInc2 .- centerEdge;

        #compare directions and switch inc if needed
        if dot(currentNormal, globalNormal)<0
            inc=inc[[2,1]];
            switched = true
        end
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #t1(v::Array{Float64,1})=transformation(m,coord1,v[1],v[2]);
        #t2(v::Array{Float64,1})=transformation(m,coord2,v[1],v[2]);

        if switched
            coordve= @views m.geometry.coordinates[:,incv[offv[e]:(offv[e]+1)]];
            if h1 #Fall periodische Randkante
                coordv= @views m.geometry.coordinates[:,incv[offv[e1]:(offv[e1]+1)]];
            else #Fall innere Kante
                coordv=coordve
            end
        else
            coordv= @views m.geometry.coordinates[:,incv[offv[e]:(offv[e]+1)]];
            if h1 #Fall periodische Randkante
                coordve= @views m.geometry.coordinates[:,incv[offv[e1]:(offv[e1]+1)]];
            else #Fall innere Kante
                coordve=coordv
            end
        end
        eT1=Array{Float64,1}();
        eT2=Array{Float64,1}();
        #coordvn1=Array{Float64,2}(undef,m.geometry.dim,mt);
        #coordvn2=Array{Float64,2}(undef,m.geometry.dim,mt);
        #for i in 1:mt
        #    coordvn1[:,i]=t1(coordref[:,i])
        #    coordvn2[:,i]=t2(coordref[:,i])
        #end
        coord1= m.geometry.coordinates[:,incf[offf[inc[1]]:(offf[inc[1]+1]-1)]];
        coord2= m.geometry.coordinates[:,incf[offf[inc[2]]:(offf[inc[2]+1]-1)]];

        v1=findall(coordve[:,1],coord1,1e-10);
        sort!(append!(v1, findall(coordve[:,2],coord1,1e-10)))
        eT1=edgeTypes[v1]
        v2=findall(coordv[:,1],coord2,1e-10);
        sort!(append!(v2, findall(coordv[:,2],coord2,1e-10)))
        eT2=edgeTypes[v2]
        globalNumVf=l2g(degFVf,inc[1])
        rb=refBound[v1]
        for j in 1:length(rb)
            if rb[j]==1.0
                push!(globv,globalNumVf[j]);
                zo+=1;
            end
        end
        push!(off,zo);
        append!(cells,inc)
        push!(edgeType,eT1)
        push!(edgeType,eT2)
        push!(edgeNum,e)
        push!(p.data.compoundData.isEdgePeriodic,h1)
    end
    p.data.edgeData=[edgeNum,cells,edgeType,globv,off];
    return nothing;
end

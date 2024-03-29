function setEdgeData!(p::femProblem, compVf::Symbol)
    m=p.mesh;
    degFVf=p.degFBoundary[p.femType[compVf][1]];
    mt=m.meshType;
    refBound=getElementProperties(mt,p.femType[compVf][1],m.geometry.dim);
    if mt==4
        edgeTypes=Dict([1,2]=>1,[2,3]=>2,[3,4]=>3,[1,4]=>4)
        coordref=[0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0]
    else
        edgeTypes=Dict([1,2]=>1,[2,3]=>2,[1,3]=>3)
        coordref=[0.0 1.0 0.0; 0.0 0.0 1.0]
    end

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
    zo=1;
    for e in 1:m.topology.size[m.topology.dim]
        off1=offe[e];
        off2=offe[e+1];
        h1=false;
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

        coord1= @views m.geometry.coordinates[:,incf[offf[inc[1]]:(offf[inc[1]]+mt-1)]];
        coord2= @views m.geometry.coordinates[:,incf[offf[inc[2]]:(offf[inc[2]]+mt-1)]];

        t1(v::Array{Float64,1})=transformation(m,coord1,v[1],v[2]);
        t2(v::Array{Float64,1})=transformation(m,coord2,v[1],v[2]);


        coordv= @views m.geometry.coordinates[:,incv[offv[e]:(offv[e]+1)]];
        if h1 #Fall periodische Randkante
            coordve= @views m.geometry.coordinates[:,incv[offv[e1]:(offv[e1]+1)]];
        else #Fall innere Kante
            coordve=coordv
        end

        eT1=Array{Float64,1}();
        eT2=Array{Float64,1}();
        coordvn1=Array{Float64,2}(undef,m.geometry.dim,mt);
        coordvn2=Array{Float64,2}(undef,m.geometry.dim,mt);
        for i in 1:mt
            coordvn1[:,i]=t1(coordref[:,i])
            coordvn2[:,i]=t2(coordref[:,i])
        end
        atol=1e-5
        v1=findall(coordve[:,1],coordvn1,atol);
        sort!(append!(v1, findall(coordve[:,2],coordvn1,atol)))
        eT1=edgeTypes[v1]
        v2=findall(coordv[:,1],coordvn2,atol);
        sort!(append!(v2, findall(coordv[:,2],coordvn2,atol)))
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
    end

    p.data.edgeData=[edgeNum,cells,edgeType,globv,off];
    return nothing;
end

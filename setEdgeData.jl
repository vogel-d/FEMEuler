function setEdgeData!(p::femProblem, compVf::Symbol)
    m=p.mesh;
    equals=p.equals;
    degFVf=p.degFBoundary[p.femType[compVf][1]];
    mt=m.meshType;
    mt==4 ?  normal=[0.5 0.0 1.0 1.0 0.5 0.0 0.0 -1.0; 0.0 -1.0 0.5 0.0 1.0 1.0 0.5 0.0] : normal=[0.5 0.0 0.5 1/sqrt(2) 0.0 -1.0; 0.0 -1.0 0.5 1/sqrt(2) 0.5 0.0];
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
    for e in 1:m.topology.size[m.topology.D]
        off1=offe[e];
        off2=offe[e+1];
        h1=false;
        if (off2-off1)==1
            e1=equals[e];
            if e1!=0
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
        mv=(1/size(coordv,2)).*[sum(coordv[1,:]), sum(coordv[2,:])];
        n1=Array{Float64,1}();
        n2=Array{Float64,1}();

        if h1 #Fall periodische Randkante
            coordve= @views m.geometry.coordinates[:,incv[offv[e1]:(offv[e1]+1)]];
            mva=(1/size(coordve,2)).*[sum(coordve[1,:]), sum(coordve[2,:])];
            mvb=mv;
        else #Fall innere Kante
            mva=mvb=mvc=mv;
        end

        for i in 1:mt
            m1=t1(normal[:,2*i-1]);
            m2=t2(normal[:,2*i-1]);
            if isapprox(m1,mva)
                n1=normal[:,2*i];
            end

            if isapprox(m2,mvb)
                n2=normal[:,2*i];
            end

        end
        globalNumVf=l2g(degFVf,inc[1])
        for i in 1:mt
            rb=t1(degFVf.referenceBoundary[i,1:2]);
            if isapprox(rb,mva)
                d=degFVf.referenceBoundary[i,3:end];
                for j in 1:length(d)
                    if d[j]==1.0
                        push!(globv,globalNumVf[j]);
                        zo+=1;
                    end
                end
                break;
            end
        end
        push!(off,zo);
        append!(cells,inc)
        push!(edgeType,getEdgeType(n1))
        push!(edgeType,getEdgeType(n2))
        push!(edgeNum,e)
    end
    p.edgeData=[edgeNum,cells,edgeType,globv,off];
    return nothing;
end

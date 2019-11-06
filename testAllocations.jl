
m=generateRectMesh(100,100);
#meshConnectivity!(m,0,2);
#incidence=m.topology.incidence["02"]
#offset=m.topology.offset["02"]

val=rand(m.topology.size[1]);
#phi=Array{Array{Float64,2},1}(undef,4);

kubPoints, kubWeights=getKub(g, m.meshType);
boundary=getBoundary(m);
b=Set{Int64}(collect(keys(boundary)));
degF=degF(m,:P1,boundary,b,kubPoints);
phi=degF.phi;

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,1},1},m::mesh)
    for k in 1:m.topology.size[3]
        ind=getIndices(incidence,offset,k);
        w=zeros(sk);
        for i in 1:length(ind)
            w+=val[ind[i]]* phi[i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,1},1})
    w=zeros(sk);
    for k in 1:m.topology.size[3]
        ind=getIndices(incidence,offset,k);
        fill!(w,0.0);
        for i in 1:length(ind)
            w+=wval[ind[i]]*phi[i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,1},1})
    w=zeros(sk);
    for k in 1:m.topology.size[3]
        ind=getIndices(incidence,offset,k);
        fill!(w,0.0);
        for i in 1:length(ind)
            @. w+= wval[ind[i]]* phi[i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,1},1})
    w=zeros(sk);
    for k in 1:m.topology.size[3]
        getIndices!(ind,incidence,offset,k);
        fill!(w,0.0);
        for i in 1:length(ind)
            @. w+= wval[ind[i]]* phi[i];
        end
    end
end


getIndices(incidence::Array{Int64,1},offset::Array{Int64,1},k::Int64)
    return incidence[offset[k]:offset[k+1]-1];
end

getIndices!(ind::Array{Int64,1},incidence::Array{Int64,1},offset::Array{Int64,1},k::Int64)
    ind=incidence[offset[k]:offset[k+1]-1];
    return nothing;
end
getIndices!(ind::Array{Int64,1},incidence::Array{Int64,1},offset::Array{Int64,1},k::Int64)
    ind[:]=incidence[offset[k]:offset[k+1]-1];
    return nothing;
end
getIndices!(ind::Array{Int64,1},incidence::Array{Int64,1},offset::Array{Int64,1},k::Int64)
    z=1;
    for i=offset[k]:offset[k+1]-1
        ind[z]=incidence[i];
    end
    return nothing;
end

function testAllocations(val::Array{Float64,1},phi::Array{Float64,3})
    w=zeros(sk);
    for k in 1:m.topology.size[3]
        getIndices!(ind,incidence,offset,k);
        fill!(w,0.0);
        for i in 1:length(ind)
            @. w+= wval[ind[i]]* phi[:,:,i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Float64,3})
    w=zeros(sk);
    for k in 1:m.topology.size[3]
        getIndices!(ind,incidence,offset,k);
        fill!(w,0.0);
        for i in 1:length(ind)
            for j in 1:sk
                for k in 1:sk
                    w[k,j]+= wval[ind[i]]* phi[k,j,i];
                end
            end

        end
    end
end


m=generateRectMesh(100,100);
#meshConnectivity!(m,0,2);
#incidence=m.topology.incidence["02"]
#offset=m.topology.offset["02"]

val=rand(m.topology.size[1]);
#phi=Array{Array{Float64,2},1}(undef,4);

kubPoints, kubWeights=getKub(g, m.meshType);
boundary=getBoundary(m);
b=Set{Int64}(collect(keys(boundary)));
degF=degF(m,:RT1,boundary,b,kubPoints);
phi=degF.phi;

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,1},1})
    for k in 1:1000
        ind=getIndices(incidence,offset,k);
        w1=zeros(sk);
        w2=zeros(sk);
        for i in 1:length(ind)
            w1+=val[globalNumW[i]]* phi[1,i];
            w2+=val[globalNumW[i]]*phi[2,i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,2},1})
    w1=zeros(sk);
    w2=zeros(sk);
    for k in 1:1000
        ind=getIndices(incidence,offset,k);
        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            w1+= wval[globalNumW[i]]* phi[1,i];
            w2+=wval[globalNumW[i]]*phi[2,i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,2},1})
    w1=zeros(sk);
    w2=zeros(sk);
    for k in 1:1000
        ind=getIndices(incidence,offset,k);
        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            @. w1+= wval[globalNumW[i]]* phi[1,i];
            @. w2+=wval[globalNumW[i]]*phi[2,i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Array{Float64,2},1})
    w1=zeros(sk);
    w2=zeros(sk);
    for k in 1:1000
        getIndices!(ind,incidence,offset,k);
        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            @. w1+= wval[globalNumW[i]]* phi[1,i];
            @. w2+=wval[globalNumW[i]]*phi[2,i];
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

function testAllocations(val::Array{Float64,1},phi::Array{Float64,4})
    w1=zeros(sk);
    w2=zeros(sk);
    for k in 1:1000
        getIndices!(ind,incidence,offset,k);
        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            @. w1+= wval[globalNumW[i]]* phi[1,:,:,i];
            @. w2+=wval[globalNumW[i]]*phi[2,:,:,i];
        end
    end
end

function testAllocations(val::Array{Float64,1},phi::Array{Float64,4})
    w1=zeros(sk);
    w2=zeros(sk);
    for k in 1:1000
        getIndices!(ind,incidence,offset,k);
        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            for j in 1:sk
                for k in 1:sk
                    w1[k,j]+= wval[globalNumW[i]]* phi[1,k,j,i];
                    w2[k,j]+=wval[globalNumW[i]]*phi[2,k,j,i];
                end
            end

        end
    end
end

function assembMassRho(degF::degF{3}, degFRho::degF{3}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,3);

    J=initPhi((2,2),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,2,m.meshType);

    globalNum=Array{Int64,1}(undef,size(phi,3));
    globalNumRho=Array{Int64,1}(undef,size(phiRho,3));

    cRho=zeros(sk);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);
        l2g!(globalNum,degF,k);
        l2g!(globalNumRho,degFRho,k);

        fill!(cRho,0.0);
        for i in 1:length(globalNumRho)
            @views @. cRho+=valRho[globalNumRho[i]]*phiRho[:,:,i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*cRho[l,r]*phi[l,r,i]*phi[l,r,j]*dJ[l,r];
                    end
                end
                if !isequal(currentval,0.0)
                    push!(rows,globalNum[i]);
                    push!(cols,globalNum[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

function assembMassRho(degF::degF{4}, degFRho::degF{3}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,4);

    J=initPhi((2,2),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=Array{Float64,4}(undef,size(phi,1),sk[1],sk[2],size(phi,4));
    coord=Array{Float64,2}(undef,2,m.meshType);

    globalNum=Array{Int64,1}(undef,size(phi,4));
    globalNumRho=Array{Int64,1}(undef,size(phiRho,3));

    cRho=zeros(sk);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints, phi,coord);
        l2g!(globalNum,degF,k);
        l2g!(globalNumRho,degFRho,k);

        fill!(cRho,0.0);
        for i in 1:length(globalNumRho)
            @views @. cRho+=valRho[globalNumRho[i]]*phiRho[:,:,i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*ddJ[l,r]*cRho[l,r]*(jphi[1,l,r,i]*jphi[1,l,r,j]+jphi[2,l,r,i]*jphi[2,l,r,j]);
                    end
                end
                if !isequal(currentval,0.0)
                    push!(rows,globalNum[i]);
                    push!(cols,globalNum[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

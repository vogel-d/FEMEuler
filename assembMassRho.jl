function assembMassRho(degF::degF{1}, degFRho::degF{1}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=length(phi);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    J=Array{Array{Float64,2},2}(undef,2,2);
    dJ=Array{Float64,2}(undef,sk);
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints);
        globalNum=@views l2g(degF,k);
        globalNumRho=@views l2g(degFRho,k);

        cRho=zeros(sk);
        for i in 1:length(globalNumRho)
            cRho+=valRho[globalNumRho[i]]*phiRho[i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*cRho[l,r]*phi[i][l,r]*phi[j][l,r]*dJ[l,r];
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

function assembMassRho(degF::degF{2}, degFRho::degF{1}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,2);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    J=Array{Array{Float64,2},2}(undef,2,2);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initPhi(size(phi),sk);
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints, phi);
        globalNum=@views l2g(degF,k);
        globalNumRho=@views l2g(degFRho,k);
        cRho=zeros(sk);
        for i in 1:length(globalNumRho)
            cRho+=valRho[globalNumRho[i]]*phiRho[i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*ddJ[l,r]*cRho[l,r]*(jphi[1,i][l,r]*jphi[1,j][l,r]+jphi[2,i][l,r]*jphi[2,j][l,r]);
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

function assembMassRho(degF::degF{1}, degFRho::degF{1}, valRho::Array{AbstractFloat,1}, m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=length(phi);
    n=degF.numB;

    J=initJacobi((2,2),sk);
    dJ=Array{AbstractFloat,2}(undef,sk);
    coord=Array{AbstractFloat,2}(undef,2,m.meshType);

    globalNum=Array{Int,1}(undef,length(phi));
    globalNumRho=Array{Int,1}(undef,length(phiRho));

    cRho=zeros(sk);

    rows=Int[];
    cols=Int[];
    vals=AbstractFloat[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);
        l2g!(globalNum,degF,k);
        l2g!(globalNumRho,degFRho,k);

        fill!(cRho,0.0);
        for i in 1:length(globalNumRho)
            @. cRho+=valRho[globalNumRho[i]]*phiRho[i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*cRho[l,r]*phi[i][l,r]*phi[j][l,r]*dJ[l,r];
                    end
                end
                if !isequal(currentval,0.0) || (globalNum[i]==n && globalNum[j]==n)
                    push!(rows,globalNum[i]);
                    push!(cols,globalNum[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

function assembMassRho(degF::degF{2}, degFRho::degF{1}, valRho::Array{AbstractFloat,1}, m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,2);
    n=degF.numB;

    J=initJacobi((2,2),sk);
    ddJ=Array{AbstractFloat,2}(undef,sk);
    jphi=initJacobi(size(phi),sk);
    coord=Array{AbstractFloat,2}(undef,2,m.meshType);

    globalNum=Array{Int,1}(undef,size(phi,2));
    globalNumRho=Array{Int,1}(undef,length(phiRho));

    cRho=zeros(sk);

    rows=Int[];
    cols=Int[];
    vals=AbstractFloat[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints, phi,coord);
        l2g!(globalNum,degF,k);
        l2g!(globalNumRho,degFRho,k);

        fill!(cRho,0.0);
        for i in 1:length(globalNumRho)
            @. cRho+=valRho[globalNumRho[i]]*phiRho[i];
        end

        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*ddJ[l,r]*cRho[l,r]*(jphi[1,i][l,r]*jphi[1,j][l,r]+jphi[2,i][l,r]*jphi[2,j][l,r]);
                    end
                end
                if !isequal(currentval,0.0) || (globalNum[i]==n && globalNum[j]==n)
                    push!(rows,globalNum[i]);
                    push!(cols,globalNum[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

function assembMassRho(degF::degF{1,:H1}, degFRho::degF{1,:H1}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=length(phi);
    n=degF.numB;

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNum=Array{Int64,1}(undef,length(phi));
    globalNumRho=Array{Int64,1}(undef,length(phiRho));

    cRho=zeros(sk);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.dim+1]
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
                        currentval+=kubWeights[l,r]*cRho[l,r]*phi[i][l,r]*phi[j][l,r]*abs(dJ[l,r]);
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

function assembMassRho(degF::degF{2,:H1div}, degFRho::degF{1,:H1}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,2);
    n=degF.numB;

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphi=initJacobi((m.geometry.dim, iter),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNum=Array{Int64,1}(undef,size(phi,2));
    globalNumRho=Array{Int64,1}(undef,length(phiRho));

    cRho=zeros(sk);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.dim+1]
        jacobi!(J,ddJ,jphi,m,k,kubPoints,phi,coord);
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
                        jphiDot=0.0;
                        for d in 1:m.geometry.dim
                            jphiDot+=jphi[d,i][l,r]*jphi[d,j][l,r];
                        end
                        currentval+=kubWeights[l,r]*abs(ddJ[l,r])*cRho[l,r]*jphiDot;
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

function assembMassRho(degF::degF{2,:H1xH1}, degFRho::degF{1,:H1}, valRho::Array{Float64,1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phi=@views degF.phi;
    phiRho=@views degFRho.phi;
    sk=size(kubWeights);
    iter=size(phi,2);
    n=degF.numB;

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNum=Array{Int64,1}(undef,size(phi,2));
    globalNumRho=Array{Int64,1}(undef,length(phiRho));

    cRho=zeros(sk);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.dim+1]
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
                        phiDot=0.0;
                        for d in 1:m.geometry.dim
                            phiDot+=phi[d,i][l,r]*phi[d,j][l,r];
                        end
                        currentval+=kubWeights[l,r]*abs(dJ[l,r])*cRho[l,r]*phiDot;
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

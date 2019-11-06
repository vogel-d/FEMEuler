#Interface zum Assemblieren der globalen Massenmatrix
function assembMass!(p::femProblem)
    degF=p.degFBoundary;
    s=Array{Symbol,1}();
    if p.taskRecovery
        for i in collect(keys(p.femType))
            append!(s, p.femType[i])
        end
    else
        for i in collect(keys(p.femType))
            push!(s, p.femType[i][1]);
        end
    end
    comp=Set{Symbol}(s);
    for i in comp
        massM=assembMass(p.degFBoundary[i],p.mesh,p.kubPoints,p.kubWeights);
        n=p.degFBoundary[i].num
        p.massM[i]=lu(massM[1:n,1:n]);
        p.massMBoundary[i]=lu(massM);
    end
end


#Funktion zum Assemblieren der globalen Massematrix für einen Finite-Elemente-Raum mit skalaren Ansatzfunktionen
function assembMass(degF::degF{3}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiRef=degF.phi;
    iter=size(phiRef,3);
    sk=size(kubWeights)
    J=initPhi((2,2),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,2,m.meshType);
    globalNum=Array{Int64,1}(undef,size(phiRef,3));

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints,coord);
        l2g!(globalNum,degF,k);
        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*phiRef[l,r,i]*phiRef[l,r,j]*dJ[l,r];
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

#Funktion zum Assemblieren der globalen Massematrix für einen Finite-Elemente-Raum mit vektoriellen Ansatzfunktionen
function assembMass(degF::degF{4}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiRef=degF.phi;
    iter=size(phiRef,4);
    sk=size(kubWeights)
    J=initPhi((2,2),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiRef=Array{Float64,4}(undef,size(phiRef,1),sk[1],sk[2],size(phiRef,4));
    coord=Array{Float64,2}(undef,2,m.meshType);
    globalNum=Array{Int64,1}(undef,size(phiRef,4));

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphiRef,m,k,kubPoints, phiRef,coord);
        l2g!(globalNum,degF,k);
        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        currentval+=kubWeights[l,r]*ddJ[l,r]*(jphiRef[1,l,r,i]*jphiRef[1,l,r,j]+jphiRef[2,l,r,i]*jphiRef[2,l,r,j]);
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

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
function assembMass(degF::degF{1}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    phiRef=degF.phi;
    iter=length(phiRef);
    J=Array{Array{Float64,2},2}(undef,2,2);
    dJ=Array{Float64,2}(undef,size(kubWeights));
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,dJ,m,k,kubPoints);
        gvertices=l2g(degF,k);
        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        currentval+=kubWeights[l,r]*phiRef[i][l,r]*phiRef[j][l,r]*dJ[l,r];
                    end
                end
                if !isequal(currentval,0.0)
                    push!(rows,gvertices[i]);
                    push!(cols,gvertices[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

#Funktion zum Assemblieren der globalen Massematrix für einen Finite-Elemente-Raum mit vektoriellen Ansatzfunktionen
function assembMass(degF::degF{2}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    rows=Int64[];
    cols=Int64[];
    vals=Float64[];
    phiRef=degF.phi;
    iter=size(phiRef,2);
    J=Array{Array{Float64,2},2}(undef,2,2);
    ddJ=Array{Float64,2}(undef,size(kubWeights));
    jphiRef=initPhi(size(phiRef),size(kubWeights));
    for k in 1:m.topology.size[m.topology.D+1]
        jacobi!(J,ddJ,jphiRef,m,k,kubPoints, phiRef);
        gvertices=l2g(degF,k);
        for j in 1:iter
            for i in 1:iter
                currentval=0.0;
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        currentval+=kubWeights[l,r]*ddJ[l,r]*(jphiRef[1,i][l,r]*jphiRef[1,j][l,r]+jphiRef[2,i][l,r]*jphiRef[2,j][l,r]);
                    end
                end
                if !isequal(currentval,0.0)
                    push!(rows,gvertices[i]);
                    push!(cols,gvertices[j]);
                    push!(vals,currentval);
                end
            end
        end
    end
    return sparse(rows,cols,vals);
end

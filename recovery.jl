function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1}, stencil)
    #cS=projectRecovery(p.degFBoundary[comp[1]],p.degFBoundary[comp[1]],cval,p.mesh,p.kubPoints,p.kubWeights)
    massM=assembMass(p.degFBoundary[comp[1]], p.mesh, p.kubPoints, p.kubWeights);
    cS=massM*cval;
    cR=recovery(p.degFBoundary[comp[1]],p.degFBoundary[comp[3]],cS,stencil,p.mesh,p.kubPoints,p.kubWeights)
    
    return cR
end

function recovery(degFT::degF{1,:H1}, degFR::degF{1,:H1}, cval::Array{Float64,1},
                  stencil, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=@views degFR.phi;
    nT=length(phiT)
    nR=length(phiR)
    sk=size(kubWeights);

    nS=length(stencil[1])

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nS*nT);
    globalNumR=Array{Int64,1}(undef,nR);

    lM=zeros(nS*nT,nR)
    cR=zeros(degFR.numB)

    for f in 1:m.topology.size[3]
        fill!(lM, 0.0);
        jacobi!(J,dJ,m,f,kubPoints,coord);
        for j in 1:nR
            z=1;
            for k in stencil[f]

                for i in 1:nT
                    currentval=0.0;
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            currentval+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]*phiR[j][l,r];
                        end
                    end
                    lM[z,j]+=currentval;
                    z+=1;
                end
            end
        end

        l2g!(globalNumR,degFR,f);
        l2g!(globalNumT,degFT,stencil[f]);

        #qr!(lM)
        #cR[globalNumR]=lM\cval[globalNumT];

        #qrM=qr(lM)
        cR[globalNumR]=lM\cval[globalNumT];

    end
    return cR;
end

function recovery(degFT::degF{2,:H1div}, degFR::degF{2,:H1xH1}, cval::Array{Float64,1},
                  stencil, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=@views degFR.phi;
    nT=size(phiT,2)
    nR=size(phiR,2)
    sk=size(kubWeights);

    nS=length(stencil[1])

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT=Array{Int64,1}(undef,nS*nT);
    globalNumR=Array{Int64,1}(undef,nR);

    lM=zeros(nS*nT,nR)
    cR=zeros(degFR.numB)

    for f in 1:m.topology.size[3]
        fill!(lM, 0.0);
        jacobi!(J,ddJ,jphiT,m,f,kubPoints,phiT,coord);
        for j in 1:nR
            z=1;
            for k in stencil[f]

                for i in 1:nT
                    currentval=0.0;
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*phiR[d,j][l,r];
                            end
                            currentval+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*vecdot;
                        end
                    end
                    lM[z,j]+=currentval;
                    z+=1;
                end
            end
        end

        l2g!(globalNumR,degFR,f);
        l2g!(globalNumT,degFT,stencil[f]);

        #qr!(lM)
        #cR[globalNumR]=lM\cval[globalNumT];

        #qrM=qr(lM)
        cR[globalNumR]=lM\cval[globalNumT];

    end
    return cR;
end



function recovery(p::femProblem, dim::Int64, comp::Array{Symbol,1}, cval::Array{Float64,1}) #,m::mesh,kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    m=p.mesh;
    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    degF=p.degFBoundary;
    cH=projectRecovery(degF[comp[2]],degF[comp[1]],cval,p.massMBoundary[comp[2]],m,kubPoints,kubWeights);
    cHP=projectRecovery(degF[comp[4]],degF[comp[2]],cH,p.massMBoundary[comp[4]],m,kubPoints,kubWeights);
    #=
    n=m.topology.size[3]
    dim=Val{dim}();
    cR=embed(p,comp[2],comp[3],cH,n,dim);
    cEmbed=embed(p,comp[1],comp[3],cval,n,dim);
    cHPEmbed=embed(p,comp[4],comp[3],cHP,n,dim);
    =#
    cR=projectRecovery(degF[comp[3]],degF[comp[2]],cH,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cEmbed=projectRecovery(degF[comp[3]],degF[comp[1]],cval,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);
    cHPEmbed=projectRecovery(degF[comp[3]],degF[comp[4]],cHP,p.massMBoundary[comp[3]],m,kubPoints,kubWeights);

    return cR+(cEmbed-cHPEmbed);
end

function recovery(p::femProblem, dim::Int64, comp::Array{Symbol,1}, cval::Array{Float64,1}, bcomp::Symbol)
    m=p.mesh;
    kubPoints=p.kubPoints;
    kubWeights=p.kubWeights;
    degF=p.degFBoundary;

    cH=projectRecovery(degF[comp[2]],degF[comp[1]],cval,p.massMBoundary[comp[2]],m,kubPoints,kubWeights);
    if haskey(p.boundaryValues,(bcomp,comp[2]))
        cH[degF[comp[2]].num+1:end]=p.boundaryValues[(bcomp,comp[2])];
    end
    cHP=projectRecovery(degF[comp[4]],degF[comp[2]],cH,p.massMBoundary[comp[4]],m,kubPoints,kubWeights);
    if haskey(p.boundaryValues,(bcomp,comp[4]))
        cH[degF[comp[4]].num+1:end]=p.boundaryValues[(bcomp,comp[4])];
    end

    n=m.topology.size[3]

    dim=Val{dim}();
    cR=embed(p,comp[2],comp[3],cH,n,dim);
    cEmbed=embed(p,comp[1],comp[3],cval,n,dim);
    cHPEmbed=embed(p,comp[4],comp[3],cHP,n,dim);

    return cR+(cEmbed-cHPEmbed);
end

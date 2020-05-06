function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1})
    cS=assembRecovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],cval,p.mesh,p.kubPoints,p.kubWeights)
    cR=recovery(p.degFBoundary[comp[2]],comp[3],cS,p.stencil,p.stencilBoundary,p.mesh,p.kubPoints,p.kubWeights)
    return cR
end

function recovery(degFT::degF{1,:H1}, recoverySpace::Symbol, cval::Array{Float64,1},
                  stencil::Array{Array{Int,1},1}, stencilBoundary::SparseMatrixCSC{Int,Int}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    #phiR=@views degFR.phi;

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    phiR=getPhiRecovery([0.0,0.0],Val(recoverySpace));
    nR=length(phiR)
    nT=length(phiT)
    #nR=div(length(phiR),nf)

    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    dxy= @. (m.geometry.r-m.geometry.l)/m.topology.n

    globalNumT=Array{Int64,1}(undef,nT);
    lcval=Float64[]
    cR=Float64[];
    zf=0;
    for f in 1:nf
        #jacobi!(J,dJ,m,f,kubPoints,coord);
        fcoord=@views mcoord[:,inc[off[f]:off[f+1]-1]]
        n=transformation(m,fcoord,0.5,0.5);
        t1,t2=getTangentialPlane(n)
        #phiR=getPhiRecovery(n,Val(recoverySpace));

        nS=length(stencil[f])
        lM=zeros(nS*nT,nR)
        z=1;
        for k in stencil[f]
            jacobi!(J,dJ,m,k,kubPoints,coord);
            k==f ? w=10.0^8 : w=1.0;
            b = @. (k==abs(stencilBoundary[:,f]))
            if iszero(b)
                kcoord=@views mcoord[:,inc[off[k]:off[k+1]-1]]
            else
                kcoord=similar(fcoord);
                if b[1]==1
                    if stencilBoundary[1,f]<0
                        dir=2; inds1=[4,3]; inds2=[1,2]
                    else
                        dir=2; inds1=[1,2]; inds2=[4,3]
                    end
                else
                    if stencilBoundary[2,f]<0
                        dir=1; inds1=[2,3]; inds2=[1,4]
                    else
                        dir=1; inds1=[1,4]; inds2=[2,3]
                    end
                end
                kcoord[:,inds1]=fcoord[:,inds2]
                kcoord[:,inds2]=fcoord[:,inds2]
                @. kcoord[dir,inds2]=fcoord[dir,inds2]-dxy[dir]
            end

            for i in 1:nT
                for j in 1:nR
                    currentval=0.0;
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            currentval+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]* #phiR[zf+j][l,r]
                                #phiR[j](transformation(m,kcoord,kubPoints[1,l],kubPoints[2,r]));
                                phiR[j](transformRecoveryCoord(n,t1,t2,transformation(m,kcoord,kubPoints[1,l],kubPoints[2,r])));
                        end
                    end
                    lM[z,j]+=w*currentval;
                end
                z+=1;
            end
            l2g!(globalNumT,degFT,k)
            append!(lcval,w.*cval[globalNumT])
        end

        qrM=qr(lM, Val(true))
        append!(cR,qrM\lcval);
        empty!(lcval);
        zf+=nR;
    end
    return cR;
end

function recovery(p::femProblem, comp::Array{Symbol,1}, cval::Array{Float64,1}, stencil)
    cS=assembRecovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[1]],cval,p.mesh,p.kubPoints,p.kubWeights)
    cR=recovery(p.degFBoundary[comp[2]],p.degFBoundary[comp[3]],cS,stencil,p.mesh,p.kubPoints,p.kubWeights)
    return cR
end

#=
function recovery(degFT::degF{1,:H1}, degFR::degF{1,:H1}, cval::Array{Float64,1},
                  stencil, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=@views degFR.phi;
    nT=length(phiT)
    nR=length(phiR)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumR=Array{Int64,1}(undef,nR);
    globalNumT=Array{Int64,1}(undef,nT);

    lcval=Float64[]
    cR=zeros(degFR.numB)
    for f in 1:m.topology.size[3]
        nS=length(stencil[f])
        lM=zeros(nS*nT,nR)

        jacobi!(J,dJ,m,f,kubPoints,coord);
        for j in 1:nR
            z=1;
            for k in stencil[f]
                k==f ? w=10^8 : w=1;
                for i in 1:nT
                    currentval=0.0;
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            currentval+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]*phiR[j][l,r];
                        end
                    end
                    lM[z,j]+=w*currentval;
                    z+=1;
                end
            end
        end

        l2g!(globalNumR,degFR,f);

        for c in stencil[f]
            l2g!(globalNumT,degFT,c)
            f==c ? w=10^8 : w=1
            append!(lcval,w.*cval[globalNumT])
        end

        qrM=qr(lM, Val(true))
        cR[globalNumR]=qrM\lcval;

        empty!(lcval);
    end
    return cR;
end
=#
function recovery(degFT::degF{2,:H1div}, degFR::degF{2,:H1xH1}, cval::Array{Float64,1},
                  stencil, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=@views degFR.phi;
    nT=size(phiT,2)
    nR=size(phiR,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumR=Array{Int64,1}(undef,nR);
    globalNumT=Array{Int64,1}(undef,nT);

    lcval=Float64[]
    cR=zeros(degFR.numB)
    for f in 1:m.topology.size[3]
        nS=length(stencil[f])
        lM=zeros(nS*nT,nR)
        jacobi!(J,ddJ,jphiT,m,f,kubPoints,phiT,coord);
        for j in 1:nR
            z=1;
            for k in stencil[f]
                k==f ? w=10^8 : w=1;
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
                    lM[z,j]+=w*currentval;
                    z+=1;
                end
            end
        end

        l2g!(globalNumR,degFR,f);

        for c in stencil[f]
            l2g!(globalNumT,degFT,c)
            f==c ? w=10^8 : w=1
            append!(lcval,w.*cval[globalNumT])
        end

        qrM=qr(lM, Val(true))
        cR[globalNumR]=qrM\lcval;

        empty!(lcval);
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

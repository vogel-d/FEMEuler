function recoveryMatrix!(p::femProblem)
    !p.taskRecovery && return nothing;
    for i in keys(p.femType)
        if length(p.femType[i])==3 && !haskey(p.recoveryM,(p.femType[i][2],p.femType[i][3]))
            p.recoveryM[(p.femType[i][2],p.femType[i][3])]=recoveryMatrix(p.degFBoundary[p.femType[i][2]],p.femType[i][3],p.stencil,p.stencilBoundary,p.mesh,p.kubPoints,p.kubWeights)
        end
    end
    return nothing
end

function recoveryMatrix(degFT::degF{1,:H1}, recoverySpace::Symbol, stencil::Array{Array{Int,1},1},
                  stencilBoundary::SparseMatrixCSC{Int,Int}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=getPhiRecovery(Val(recoverySpace));

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    nR=length(phiR)
    nT=length(phiT)

    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    dxy= @. (m.geometry.r-m.geometry.l)/m.topology.n
    recoveryM=Array{QRPivoted{Float64,Array{Float64,2}},1}(undef,nf)
    for f in 1:nf
        fcoord=@views mcoord[:,inc[off[f]:off[f+1]-1]]
        n=transformation(m,fcoord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

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
                            currentval+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]*
                                phiR[j](transformRecoveryCoord(n,t1,t2,transformation(m,kcoord,kubPoints[1,l],kubPoints[2,r])));
                        end
                    end
                    lM[z,j]+=w*currentval;
                end
                z+=1;
            end
        end
        recoveryM[f]=qr(lM, Val(true))
    end
    return recoveryM;
end

function recoveryMatrix(degFT::degF{2,:H1div}, recoverySpace::Symbol,
                  stencil::Array{Array{Int,1},1}, stencilBoundary::SparseMatrixCSC{Int,Int}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    phiR=getPhiRecovery(Val(recoverySpace));

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    nT=size(phiT,2)
    nR=size(phiR,2)
    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    dxy= @. (m.geometry.r-m.geometry.l)/m.topology.n

    globalNumT=Array{Int64,1}(undef,nT);
    lcval=Float64[]
    cR=Float64[];
    for f in 1:nf
        fcoord=@views mcoord[:,inc[off[f]:off[f+1]-1]]
        n=transformation(m,fcoord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        nS=length(stencil[f])
        lM=zeros(nS*nT,nR)
        z=1;
        for k in stencil[f]
            jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,coord);
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
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*phiR[d,j](transformRecoveryCoord(n,t1,t2,transformation(m,kcoord,kubPoints[1,l],kubPoints[2,r])));
                            end
                            currentval+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*vecdot;
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
    end
    return cR;
end

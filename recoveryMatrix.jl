function recoveryMatrix!(p::femProblem)
    !p.taskRecovery && return nothing;
    for i in keys(p.femType)
        if length(p.femType[i])==3 && !haskey(p.recoveryM,(p.femType[i][2],p.femType[i][3]))
            p.recoveryM[(p.femType[i][2],p.femType[i][3])]=recoveryMatrix(p.degFBoundary[p.femType[i][2]],p.degFBoundary[p.femType[i][3]],p.femType[i][3],p.stencil,p.stencilBoundary,p.mesh,p.kubPoints,p.kubWeights)
        end
    end
    return nothing
end

function recoveryMatrix(degFT::degF{1,:H1}, degFR::degF{1,:H1}, recoverySpace::Symbol, stencil::Array{Array{Int,1},1},
                  stencilBoundary::SparseMatrixCSC{Int,Int}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})
    phiT=@views degFT.phi;
    nR=getPhiRecoveryLength(recoverySpace)
    nT=length(phiT)
    phiR=zeros(nR);

    globalNumT=zeros(Int,nT);
    globalNumR=zeros(Int,nR);

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    kubP=zeros(m.geometry.dim)

    ind=[1,2,3]
    Ji=zeros(3,3)
    x=zeros(3)
    kpoint=zeros(m.geometry.dim)

    dxy= @. (m.geometry.r-m.geometry.l)/m.topology.n

    lM=zeros(nT,nR)
    recoveryM=Array{Any,1}(undef,nf)
    #recoveryM=Array{LinearAlgebra.QRCompactWY{Float64,Array{Float64,2}},1}(undef,nf)
    for f in 1:nf
        fcoord=@views mcoord[:,inc[off[f]:off[f+1]-1]]

        nS=length(stencil[f])
        rows=zeros(Int,nT*nR*nS)
        cols=zeros(Int,nT*nR*nS)
        vals=zeros(nT*nR*nS)
        l2g!(globalNumR,degFR,f);
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
            l2g!(globalNumT,degFT,k);
            fill!(lM,0.0)
            for r in 1:sk[2]
                for l in 1:sk[1]
                    transformation!(kubP,m,kcoord,kubPoints[1,l],kubPoints[2,r])
                    ksi,eta=intersect(fcoord,kubP,Ji,x);
                    getPhiRecovery!(phiR,[ksi,eta])
                    for i in 1:nT
                        for j in 1:nR
                            lM[i,j]+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]*phiR[j];
                        end
                    end
                end
            end
            for i=1:nT,j=1:nR
                vals[z]=w*lM[i,j]
                rows[z]=globalNumT[i]
                cols[z]=globalNumR[j]
                z+=1
            end
            #=
            for i in 1:nT
                for j in 1:nR
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            transformation!(kubP,m,kcoord,kubPoints[1,l],kubPoints[2,r])
                            ksi,eta=intersect(fcoord,kubP,Ji,x);
                            getPhiRecovery!(phiR,[ksi,eta])
                            lM[i,j]+=kubWeights[l,r]*abs(dJ[l,r])*phiT[i][l,r]*phiR[j];
                        end
                    end
                    vals[z]=w*lM[i,j]
                    rows[z]=globalNumT[i]
                    cols[z]=globalNumR[j]
                    z+=1
                end
            end
            =#
        end
        M=sparse(rows,cols,vals)
        unique!(rows)
        unique!(cols)
        l2g!(globalNumT,degFT,f);
        for i in 1:length(globalNumT)
            for j in 1:length(rows)
                if globalNumT[i]==rows[j]
                    rows[j]=rows[i]
                    rows[i]=globalNumT[i]
                    break
                end
            end
        end
        recoveryM[f]=(qr(Matrix(M[rows,cols]), Val(true)),rows)
    end
    return recoveryM;
end

function recoveryMatrix(degFT::degF{2,:H1div}, degFR::degF{2,:H1xH1}, recoverySpace::Symbol,
                  stencil::Array{Array{Int,1},1}, stencilBoundary::SparseMatrixCSC{Int,Int}, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})

    phiT=@views degFT.phi;
    nR=getPhiRecoveryLength(recoverySpace)
    nT=size(phiT,2)
    phiR=zeros(m.geometry.dim,nR);

    globalNumT=zeros(Int,nT);
    globalNumR=zeros(Int,nR);

    nf=m.topology.size[3]
    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];

    sk=size(kubWeights);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,nT),sk);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    n=zeros(m.geometry.dim)
    kubP=zeros(m.geometry.dim)

    ind=[1,2,3]
    Ji=zeros(3,3)
    x=zeros(3)
    kpoint=zeros(m.geometry.dim)

    dxy= @. (m.geometry.r-m.geometry.l)/m.topology.n

    lM=zeros(nT,nR)
    recoveryM=Array{Any,1}(undef,nf)
    for f in 1:nf
        fcoord=@views mcoord[:,inc[off[f]:off[f+1]-1]]

        nS=length(stencil[f])
        rows=zeros(Int,nT*nR*nS)
        cols=zeros(Int,nT*nR*nS)
        vals=zeros(nT*nR*nS)
        normB=zeros(nS,nR)
        l2g!(globalNumR,degFR,f);
        z=1;
        zk=1;
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
            l2g!(globalNumT,degFT,k);
            fill!(lM,0.0)

            for r in 1:sk[2]
                for l in 1:sk[1]
                    transformation!(kubP,m,kcoord,kubPoints[1,l],kubPoints[2,r])
                    ksi,eta=intersect(fcoord,kubP,Ji,x);
                    getPhiRecovery!(phiR,[ksi,eta])
                    for i in 1:nT
                        for j in 1:nR
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*phiR[d,j];
                            end
                            lM[i,j]+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*vecdot;
                        end
                    end
                end
            end
            for i=1:nT,j=1:nR
                vals[z]=w*lM[i,j]
                rows[z]=globalNumT[i]
                cols[z]=globalNumR[j]
                z+=1
            end
            #=
            for i in 1:nT
                for j in 1:nR
                    for r in 1:sk[2]
                        for l in 1:sk[1]
                            transformation!(kubP,m,kcoord,kubPoints[1,l],kubPoints[2,r])
                            ksi,eta=intersect(fcoord,kubP,Ji,x);
                            getPhiRecovery!(phiR,[ksi,eta])
                            vecdot=0.0
                            for d in 1:m.geometry.dim
                                vecdot+=jphiT[d,i][l,r]*phiR[d,j];
                            end
                            lM[i,j]+=kubWeights[l,r]*ddJ[l,r]/abs(ddJ[l,r])*vecdot;
                        end
                    end
                    vals[z]=w*lM[i,j]
                    rows[z]=globalNumT[i]
                    cols[z]=globalNumR[j]
                    z+=1
                end
            end
            =#
            if m.geometry.dim==3
                transformation!(n,m,fcoord,0.5,0.5);
                transformation!(kubP,m,kcoord,0.5,0.5)
                ksi,eta=intersect(fcoord,kubP,Ji,x);
                getPhiRecovery!(phiR,[ksi,eta])
                normB[zk,:]=n'*phiR
                zk+=1
            end
        end
        Ms=sparse(rows,cols,vals)
        unique!(rows)
        unique!(cols)
        l2g!(globalNumT,degFT,f);
        for i in 1:length(globalNumT)
            for j in 1:length(rows)
                if globalNumT[i]==rows[j]
                    rows[j]=rows[i]
                    rows[i]=globalNumT[i]
                    break
                end
            end
        end
        #=
        M=Matrix(Ms[rows,cols])
        for k in stencil[f]
            kcoord=@views mcoord[:,inc[off[k]:off[k+1]-1]]
            transformation!(n,m,fcoord,0.5,0.5);
            transformation!(kubP,m,kcoord,0.5,0.5)
            ksi,eta=intersect(fcoord,kubP,Ji,x);
            getPhiRecovery!(phiR,[ksi,eta])
            M=vcat(M,n'*phiR)
        end
        =#
        if m.geometry.dim==3
            recoveryM[f]=(qr(vcat(Matrix(Ms[rows,cols]),normB), Val(true)),rows)
        else
            recoveryM[f]=(qr(Matrix(Ms[rows,cols]), Val(true)),rows)
        end
    end
    return recoveryM;
end

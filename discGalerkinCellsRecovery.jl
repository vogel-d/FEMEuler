function discGalerkinCellsR!(M::Array{Float64,2},
                            degFT::degF{1,:H1}, gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            phiW::Array{Function,1}, wval::Array{Float64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    nW=length(phiW)

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    zf=1
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        n=transformation(m,coord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        fill!(w,0.0);
        for i in 1:nW
            for r in 1:size(kubWeights,2)
                for l in 1:size(kubWeights,1)
                    w[l,r]+=wval[zf]*phiW[i](transformRecoveryCoord(n,t1,t2,transformation(m,coord,kubPoints[1,l],kubPoints[2,r])));
                end
            end
            zf+=1
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        phiFgradphiT=0.0;
                        for d in 1:m.topology.dim
                            phiFgradphiT+=phiF[d,j][l,r]*gradphiT[d,i][l,r]
                        end
                        z+=fval[gj]*kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*w[l,r]*phiFgradphiT;
                    end
                end
            end
            M[gi]+=z;
        end
    end

    return nothing;
end

function discGalerkinCellsR!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1}, gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div}, phiF::Array{Array{Float64,2},2}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            phiW::Array{Function,1}, wval::Array{Float64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    nW=length(phiW)
    nT=length(globalNumT);
    nF=length(globalNumF);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nT,nF);
    zf=1
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        n=transformation(m,coord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        fill!(w,0.0);
        for i in 1:nW
            for r in 1:size(kubWeights,2)
                for l in 1:size(kubWeights,1)
                    w[l,r]+=wval[zf]*phiW[i](transformRecoveryCoord(n,t1,t2,transformation(m,coord,kubPoints[1,l],kubPoints[2,r])));
                end
            end
            zf+=1
        end


        fill!(lM,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        phiFgradphiT=0.0;
                        for d in 1:m.topology.dim
                            phiFgradphiT+=phiF[d,j][l,r]*gradphiT[d,i][l,r]
                        end
                        lM[i,j]+=kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*w[l,r]*phiFgradphiT;
                    end
                end
            end
        end
        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            for j in 1:length(globalNumF)
                push!(rows,gi);
                push!(cols,globalNumF[j]);
                push!(vals,lM[i,j]);
            end
        end
    end

    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            #degFW::degF{2,S} where S,phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi(size(phiT),sk);

    w1=zeros(sk);
    w2=zeros(sk);
    gradw11=zeros(sk);
    gradw12=zeros(sk);
    gradw21=zeros(sk);
    gradw22=zeros(sk);

    for k in 1:m.topology.size[3]
        jacobi!(ddJ,jphiT,m,k,kubPoints,phiT,coord);

        l2g!(globalNumW,degFW,k);

        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            @. w1+=wval[globalNumW[i]]*phiW[1,i];
            @. w2+=wval[globalNumW[i]]*phiW[2,i];
        end

        fill!(gradw11,0.0);
        fill!(gradw12,0.0);
        fill!(gradw21,0.0);
        fill!(gradw22,0.0);
        zg=0;
        for i in 1:size(phiW,2)
            @. gradw11+=wval[globalNumW[i]]*gradphiW[1,1+zg];
            @. gradw12+=wval[globalNumW[i]]*gradphiW[1,2+zg];
            @. gradw21+=wval[globalNumW[i]]*gradphiW[2,1+zg];
            @. gradw22+=wval[globalNumW[i]]*gradphiW[2,2+zg];
            zg+=2;
        end


        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        z+=fval[gj]*kubWeights[l,r]*(ddJ[l,r]^3/abs(ddJ[l,r]))*(jphiT[1,i][l,r]*(dphiF[j][l,r]*w1[l,r]+gradw11[l,r]*phiF[1,j][l,r]+gradw12[l,r]*phiF[2,j][l,r])+jphiT[2,i][l,r]*(dphiF[j][l,r]*w2[l,r]+gradw21[l,r]*phiF[1,j][l,r]+gradw22[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
            M[gi]-=z;
        end

    end
    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},phiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1xH1},phiW::Array{Array{Float64,2},2}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initJacobi((m.geometry.dim,size(phiT,2)),sk);

    w=[zeros(sk) for d in 1:m.geometry.dim]
    gradw=Array{Array{Float64,2},2}(undef,m.geometry.dim,m.topology.dim);
    for k in 1:m.topology.dim
        for j in 1:m.geometry.dim
            gradw[j,k]=zeros(sk);
        end
    end

    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,coord);
        l2g!(globalNumW,degFW,k);

        for d in 1:m.geometry.dim
            fill!(w[d], 0.0)
            for i in 1:length(globalNumW)
                @. w[d]+=wval[globalNumW[i]]*phiW[d,i];
            end
        end

        for dt in 1:m.topology.dim
            for dg in 1:m.geometry.dim
                fill!(gradw[dg,dt], 0.0)
                zg=0
                for i in 1:size(phiW,2)
                    @. gradw[dg,dt]+=wval[globalNumW[i]]*gradphiW[dg,dt+zg];
                    zg+=2;
                end
            end
        end


        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        vecdot=0.0;
                        for d in 1:m.geometry.dim
                            gradwPhiF=0.0
                            for dt in 1:m.topology.dim
                                gradwPhiF+=gradw[d,dt][l,r]*phiF[dt,j][l,r]
                            end
                            vecdot+=jphiT[d,i][l,r]*(dphiF[j][l,r]*w[d][l,r]+gradwPhiF)
                        end
                        z+=fval[gj]*kubWeights[l,r]*abs(ddJ[l,r])*vecdot
                    end
                end
            end
            M[gi]-=z;
        end

    end
    return nothing;
end

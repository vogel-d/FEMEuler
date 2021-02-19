function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{1,:H1},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        jacobi!(J,dJ,m,k,kubPoints,coord);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            @. w+=wval[globalNumW[i]]*phiW[i];
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


function discGalerkinCells!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})


    sk=size(kubWeights);

    nT=length(globalNumT);
    nF=length(globalNumF);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nT,nF);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        jacobi!(J,dJ,m,k,kubPoints,coord);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            @. w+=wval[globalNumW[i]]*phiW[i];
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

        for i in 1:nT
            gi=globalNumT[i];
            for j in 1:nF
                push!(rows,gi);
                push!(cols,globalNumF[j]);
                push!(vals,lM[i,j]);
            end
        end
    end

    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);

    w1=zeros(sk);
    w2=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiW=initJacobi((m.geometry.dim,size(phiW,2)),sk);

    w=[zeros(sk) for d in 1:m.geometry.dim]

    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiW,m,k,kubPoints,phiW,coord);

        l2g!(globalNumW,degFW,k);
        for d in 1:m.geometry.dim
            fill!(w[d], 0.0)
            for i in 1:length(globalNumW)
                @. w[d]+=wval[globalNumW[i]]*jphiW[d,i];
            end
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);
        zg=0;
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        val=0.0;
                        for d in 1:m.geometry.dim
                            jgradphiTphiF=0.0;
                            for dt in 1:m.topology.dim
                                gradphiTphiF=0.0;
                                for dti in 1:m.topology.dim
                                    gradphiTphiF+=gradphiT[dt,zg+dti][l,r]*phiF[dti,j][l,r]
                                end
                                jgradphiTphiF+=J[d,dt][l,r]*gradphiTphiF
                            end
                            val+=jgradphiTphiF*w[d][l,r]
                        end
                        z+=fval[gj]*kubWeights[l,r]*(ddJ[l,r]^3)/abs(ddJ[l,r])*val
                    end
                end
            end
            M[gi]+=z;
            zg+=2;
        end
    end

    return nothing;
end





function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{1,:H1},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1xH1},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        jacobi!(J,dJ,m,k,kubPoints,coord);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            @. w+=wval[globalNumW[i]]*phiW[i];
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


function discGalerkinCells!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})


    sk=size(kubWeights);

    nT=length(globalNumT);
    nF=length(globalNumF);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nT,nF);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        jacobi!(J,dJ,m,k,kubPoints,coord);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            @. w+=wval[globalNumW[i]]*phiW[i];
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

        for i in 1:nT
            gi=globalNumT[i];
            for j in 1:nF
                push!(rows,gi);
                push!(cols,globalNumF[j]);
                push!(vals,lM[i,j]);
            end
        end
    end

    return nothing;
end


function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2,:H1div},gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2,:H1div},phiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);

    w1=zeros(sk);
    w2=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    ddJ=Array{Float64,2}(undef,sk);
    jphiW=initJacobi((m.geometry.dim,size(phiW,2)),sk);

    w=[zeros(sk) for d in 1:m.geometry.dim]

    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiW,m,k,kubPoints,phiW,coord);

        l2g!(globalNumW,degFW,k);
        for d in 1:m.geometry.dim
            fill!(w[d], 0.0)
            for i in 1:length(globalNumW)
                @. w[d]+=wval[globalNumW[i]]*jphiW[d,i];
            end
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);
        zg=0;
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        val=0.0;
                        for d in 1:m.geometry.dim
                            jgradphiTphiF=0.0;
                            for dt in 1:m.topology.dim
                                gradphiTphiF=0.0;
                                for dti in 1:m.topology.dim
                                    gradphiTphiF+=gradphiT[dt,zg+dti][l,r]*phiF[dti,j][l,r]
                                end
                                jgradphiTphiF+=J[d,dt][l,r]*gradphiTphiF
                            end
                            val+=jgradphiTphiF*w[d][l,r]
                        end
                        z+=fval[gj]*kubWeights[l,r]*(ddJ[l,r]^3)/abs(ddJ[l,r])*val
                    end
                end
            end
            M[gi]+=z;
            zg+=2;
        end
    end

    return nothing;
end

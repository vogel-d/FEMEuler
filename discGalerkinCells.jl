function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);

    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        jacobi!(J,dJ,m,k,kubPoints,coord);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            @. w+=wval[globalNumW[i]]*phiW[i];
        end

        fill!(gradw1,0.0);
        fill!(gradw2,0.0);
        for j in 1:size(gradphiW,2)
            @. gradw1+=wval[globalNumW[j]]*gradphiW[1,j]
            @. gradw2+=wval[globalNumW[j]]*gradphiW[2,j];
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                zLoc=0.0;
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        z+=fval[gj]*kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[i][l,r]*(w[l,r]*dphiF[j][l,r]+(gradw1[l,r]*phiF[1,j][l,r]+gradw2[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
            M[gi]-=z;
        end
    end

    return nothing;
end

function discGalerkinCells!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1},phiT::Array{Array{Float64,2},1}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, dphiF::Array{Array{Float64,2},1}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{1,:H1},phiW::Array{Array{Float64,2},1}, gradphiW::Array{Array{Float64,2},2}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})


    sk=size(kubWeights);

    nT=length(phiT);
    nF=size(phiF,2);

    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);

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

        fill!(gradw1,0.0);
        fill!(gradw2,0.0);
        for j in 1:size(gradphiW,2)
            @. gradw1+=wval[globalNumW[j]]*gradphiW[1,j]
            @. gradw2+=wval[globalNumW[j]]*gradphiW[2,j];
        end

        fill!(lM,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        lM[i,j]+=kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*phiT[i][l,r]*(w[l,r]*dphiF[j][l,r]+(gradw1[l,r]*phiF[1,j][l,r]+gradw2[l,r]*phiF[2,j][l,r]));
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
                push!(vals,-lM[i,j]);
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
    jphiT=initJacobi(size(phiT),sk);
    jphiF=initJacobi(size(phiF),sk);

    w1=zeros(sk);
    w2=zeros(sk);
    gradw11=zeros(sk);
    gradw12=zeros(sk);
    gradw21=zeros(sk);
    gradw22=zeros(sk);

    for k in 1:m.topology.size[3]
        jacobi!(J,ddJ,jphiT,m,k,kubPoints,phiT,coord);
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
                        z+=fval[gj]*kubWeights[l,r]*abs(ddJ[l,r])*(jphiT[1,i][l,r]*(dphiF[j][l,r]*w1[l,r]+gradw11[l,r]*phiF[1,j][l,r]+gradw12[l,r]*phiF[2,j][l,r])+jphiT[2,i][l,r]*(dphiF[j][l,r]*w2[l,r]+gradw21[l,r]*phiF[1,j][l,r]+gradw22[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
            M[gi]-=z;
        end

    end
    return nothing;
end

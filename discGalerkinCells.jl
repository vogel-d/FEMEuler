function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{3},phiT::Array{Float64,3}, globalNumT::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, dphiF::Array{Float64,3}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{3},phiW::Array{Float64,3}, gradphiW::Array{Float64,4}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})


    sk=size(kubWeights);

    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);

    @inbounds for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            #@views @. w+=wval[globalNumW[i]]*phiW[:,:,i];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    w[l,r]+= wval[globalNumW[i]]*phiW[l,r,i];
                end
            end
        end

        fill!(gradw1,0.0);
        fill!(gradw2,0.0);
        for j in 1:length(globalNumW)
            #@views @. gradw1+=wval[globalNumW[j]]*gradphiW[1,:,:,j]
            #@views @. gradw2+=wval[globalNumW[j]]*gradphiW[2,:,:,j];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gradw1[l,r]+=wval[globalNumW[j]]*gradphiW[1,l,r,j];
                    gradw2[l,r]+=wval[globalNumW[j]]*gradphiW[2,l,r,j];
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
                        z+=fval[gj]*kubWeights[l,r]*phiT[l,r,i]*(w[l,r]*dphiF[l,r,j]+(gradw1[l,r]*phiF[1,l,r,j]+gradw2[l,r]*phiF[2,l,r,j]));
                    end
                end
            end
            M[gi]-=z;
        end
    end

    return nothing;
end

function discGalerkinCells!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{3},phiT::Array{Float64,3}, globalNumT::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, dphiF::Array{Float64,3}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            degFW::degF{3},phiW::Array{Float64,3}, gradphiW::Array{Float64,4}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2})


    sk=size(kubWeights);

    nT=size(phiT,3);
    nF=size(phiF,4);

    w=zeros(sk);
    gradw1=zeros(sk);
    gradw2=zeros(sk);

    lM=zeros(nT,nF);

    for k in 1:m.topology.size[3]
        l2g!(globalNumW,degFW,k);

        fill!(w,0.0);
        for i in 1:length(globalNumW)
            #@views @. w+=wval[globalNumW[i]]*phiW[:,:,i];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    w[l,r]+= wval[globalNumW[i]]*phiW[l,r,i];
                end
            end
        end

        fill!(gradw1,0.0);
        fill!(gradw2,0.0);
        for j in 1:length(globalNumW)
            #@views @. gradw1+=wval[globalNumW[j]]*gradphiW[1,:,:,j]
            #@views @. gradw2+=wval[globalNumW[j]]*gradphiW[2,:,:,j];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gradw1[l,r]+=wval[globalNumW[j]]*gradphiW[1,l,r,j];
                    gradw2[l,r]+=wval[globalNumW[j]]*gradphiW[2,l,r,j];
                end
            end
        end

        fill!(lM,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        lM[i,j]+=kubWeights[l,r]*phiT[l,r,i]*(w[l,r]*dphiF[l,r,j]+(gradw1[l,r]*phiF[1,l,r,j]+gradw2[l,r]*phiF[2,l,r,j]));
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
                            degFT::degF{4},phiT::Array{Float64,4}, globalNumT::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, dphiF::Array{Float64,3}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{4},phiW::Array{Float64,4}, gradphiW::Array{Float64,4}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})


    sk=size(kubWeights);

    ddJ=Array{Float64,2}(undef,sk);
    jphiT=Array{Float64,4}(undef,size(phiT,1),sk[1],sk[2],size(phiT,4));

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
            #@views @. w1+= wval[globalNumW[i]]* phiW[1,:,:,i];
            #@views @. w2+=wval[globalNumW[i]]*phiW[2,:,:,i];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    w1[l,r]+=wval[globalNumW[i]]*phiW[1,l,r,i];
                    w2[l,r]+=wval[globalNumW[i]]*phiW[2,l,r,i];
                end
            end
        end

        fill!(gradw11,0.0);
        fill!(gradw12,0.0);
        fill!(gradw21,0.0);
        fill!(gradw22,0.0);
        zg=0;
        for i in 1:length(globalNumW)
            #@views @. gradw11+=wval[globalNumW[i]]*gradphiW[1,:,:,1+zg];
            #@views @. gradw12+=wval[globalNumW[i]]*gradphiW[1,:,:,2+zg];
            #@views @. gradw21+=wval[globalNumW[i]]*gradphiW[2,:,:,1+zg];
            #@views @. gradw22+=wval[globalNumW[i]]*gradphiW[2,:,:,2+zg];
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gradw11[l,r]+=wval[globalNumW[i]]*gradphiW[1,l,r,1+zg];
                    gradw12[l,r]+=wval[globalNumW[i]]*gradphiW[1,l,r,2+zg];
                    gradw21[l,r]+=wval[globalNumW[i]]*gradphiW[2,l,r,1+zg];
                    gradw22[l,r]+=wval[globalNumW[i]]*gradphiW[2,l,r,2+zg];
                end
            end
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
                        z+=fval[gj]*kubWeights[l,r]*ddJ[l,r]^2*(jphiT[1,l,r,i]*(dphiF[l,r,j]*w1[l,r]+gradw11[l,r]*phiF[1,l,r,j]+gradw12[l,r]*phiF[2,l,r,j])+jphiT[2,l,r,i]*(dphiF[l,r,j]*w2[l,r]+gradw21[l,r]*phiF[1,l,r,j]+gradw22[l,r]*phiF[2,l,r,j]));
                    end
                end
            end
            M[gi]-=z;
        end

    end
    return nothing;
end

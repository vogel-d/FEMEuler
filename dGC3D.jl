function discGalerkinCells!(M::Array{Float64,2},
                            degFT::degF{2},phiT::Array{Float64,4}, globalNumT::Array{Int64,1},
                            degFF::degF{2},phiF::Array{Float64,4}, dphiF::Array{Float64,3}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            degFW::degF{2},phiW::Array{Float64,4}, gradphiW::Array{Float64,4}, wval::Array{Float64,1}, globalNumW::Array{Int64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})


    sk=size(kubWeights);

    ddJ=Array{Float64,2}(undef,sk);
    jphiT=Array{Float64,4}(undef,size(phiT,1),sk[1],sk[2],size(phiT,2));

    w1=zeros(sk);
    w2=zeros(sk);
    gradw11=zeros(sk);
    gradw12=zeros(sk);
    gradw21=zeros(sk);
    gradw22=zeros(sk);

    @inbounds for k in 1:m.topology.size[3]
        jacobi!(ddJ,jphiT,m,k,kubPoints,phiT,coord);
        l2g!(globalNumW,degFW,k);

        fill!(w1,0.0);
        fill!(w2,0.0);
        for i in 1:length(globalNumW)
            #@. w1+= wval[globalNumW[i]]* phiW[:,:,1,i];
            #@. w2+=wval[globalNumW[i]]*phiW[:,:,2,i];

            for l in 1:size(phiW,2)
                for k in 1:size(phiW,1)
                    w1[k,l]+= wval[globalNumW[i]]* phiW[1,k,l,i];
                    w2[k,l]+=wval[globalNumW[i]]*phiW[2,k,l,i];
                end
            end

        end

        fill!(gradw11,0.0);
        fill!(gradw12,0.0);
        fill!(gradw21,0.0);
        fill!(gradw22,0.0);
        zg=0;
        for i in 1:size(phiW,4)
            for l in 1:size(gradphiW,2)
                for k in 1:size(gradphiW,1)
                    gradw11[k,l]+=wval[globalNumW[i]]*gradphiW[1,k,l,1+zg];
                    gradw21[k,l]+=wval[globalNumW[i]]*gradphiW[2,k,l,1+zg];
                    gradw12[k,l]+=wval[globalNumW[i]]*gradphiW[1,k,l,2+zg];
                    gradw22[k,l]+=wval[globalNumW[i]]*gradphiW[2,k,l,2+zg];
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

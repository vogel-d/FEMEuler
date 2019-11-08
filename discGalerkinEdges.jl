function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{3},phiT::Array{Float64,3}, phiTtrans::Array{Float64,4}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, phiFtrans::Array{Float64,4}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{3},phiW::Array{Float64,3}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64)



    nT=size(phiT,3);
    nF=size(phiF,4);
    sk=length(quadWeights)

    w1=zeros(sk);
    w2=zeros(sk);

    lM11=zeros(nT,nF);
    lM12=zeros(nT,nF);
    lM21=zeros(nT,nF);
    lM22=zeros(nT,nF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        n=@views m.normals[:,eT1];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        phiFn1=@views phiFtrans[:,:,:,eT1];
        phiTn1=@views phiTtrans[:,:,:,eT1];
        phiWn1=@views phiWtrans[:,:,:,eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[:,:,:,eT2];
        phiTn2=@views phiTtrans[:,:,:,eT2];
        phiWn2=@views phiWtrans[:,:,:,eT2];
        kubPn2=@views nquadPoints[eT2];

        fill!(w1,0.0);
        fill!(w2,0.0);
        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        for i in 1:length(globalNumW1)
            @views @. w1+=wval[globalNumW1[i]]*phiWn1[1,:,i];
            @views @. w2+=wval[globalNumW2[i]]*phiWn2[1,:,i];
        end

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[1,r,i]*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[1,r,i]*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[1,r,i]*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[1,r,i]*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j]);
                end
            end
        end
        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];
                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end
        end
    end
    return nothing;
end

function discGalerkinEdges!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{3},phiT::Array{Float64,3}, phiTtrans::Array{Float64,4}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, phiFtrans::Array{Float64,4}, fval::Array{Float64,1}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{3},phiW::Array{Float64,3}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64)



    nT=size(phiT,3);
    nF=size(phiF,4);
    sk=length(quadWeights)

    w1=zeros(sk);
    w2=zeros(sk);

    lM11=zeros(nT,nF);
    lM12=zeros(nT,nF);
    lM21=zeros(nT,nF);
    lM22=zeros(nT,nF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        n=@views m.normals[:,eT1];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        phiFn1=@views phiFtrans[:,:,:,eT1];
        phiTn1=@views phiTtrans[:,:,:,eT1];
        phiWn1=@views phiWtrans[:,:,:,eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[:,:,:,eT2];
        phiTn2=@views phiTtrans[:,:,:,eT2];
        phiWn2=@views phiWtrans[:,:,:,eT2];
        kubPn2=@views nquadPoints[eT2];

        fill!(w1,0.0);
        fill!(w2,0.0);
        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        for i in 1:length(globalNumW1)
            @views @. w1+=wval[globalNumW1[i]]*phiWn1[1,:,i];
            @views @. w2+=wval[globalNumW2[i]]*phiWn2[1,:,i];
        end

        s=0.0;
        for i in globv
            i>length(fval) && continue;
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[1,r,i]*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[1,r,i]*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[1,r,i]*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[1,r,i]*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j]);
                end
            end
        end
        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];

                push!(rows,gi1);
                push!(cols,gj1);
                push!(vals,(+0.5-gammaLoc)*lM11[i,j]);

                push!(rows,gi1);
                push!(cols,gj2);
                push!(vals,(-0.5+gammaLoc)*lM12[i,j]);

                push!(rows,gi2);
                push!(cols,gj1);
                push!(vals,(+0.5+gammaLoc)*lM21[i,j]);

                push!(rows,gi2);
                push!(cols,gj2);
                push!(vals,(-0.5-gammaLoc)*lM22[i,j]);
            end
        end

    end
    return nothing;
end

function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{4},phiT::Array{Float64,4}, phiTtrans::Array{Float64,4}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{4},phiF::Array{Float64,4}, phiFtrans::Array{Float64,4}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{4},phiW::Array{Float64,4}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64, coord::Array{Float64,2})



    nT=size(phiT,4);
    nF=size(phiF,4);
    sk=length(quadWeights)

    J1=initPhi((2,2),sk);
    ddJ1=Array{Float64,1}(undef,sk);
    jphiWn1=Array{Float64,3}(undef,size(phiW,1),sk,size(phiW,4));
    jphiTn1=Array{Float64,3}(undef,size(phiT,1),sk,size(phiT,4));

    J2=initPhi((2,2),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiWn2=Array{Float64,3}(undef,size(phiW,1),sk,size(phiW,4));
    jphiTn2=Array{Float64,3}(undef,size(phiT,1),sk,size(phiT,4));

    w11=zeros(sk);
    w12=zeros(sk);
    w21=zeros(sk);
    w22=zeros(sk);

    lM11=zeros(nT,nF);
    lM12=zeros(nT,nF);
    lM21=zeros(nT,nF);
    lM22=zeros(nT,nF);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        n=@views m.normals[:,eT1];
        le=m.edgeLength[edgeData[1][e]];
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        phiFn1=@views phiFtrans[:,:,:,eT1];
        phiTn1=@views phiTtrans[:,:,:,eT1];
        phiWn1=@views phiWtrans[:,:,:,eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[:,:,:,eT2];
        phiTn2=@views phiTtrans[:,:,:,eT2];
        phiWn2=@views phiWtrans[:,:,:,eT2];
        kubPn2=@views nquadPoints[eT2];

        jacobi!(J1,ddJ1,jphiWn1,jphiTn1,m,inc1,kubPn1, phiWn1, phiTn1, coord);
        jacobi!(J2,ddJ2,jphiWn2,jphiTn2,m,inc2,kubPn2, phiWn2, phiTn2, coord);


        fill!(w11,0.0);
        fill!(w12,0.0);
        fill!(w21,0.0);
        fill!(w22,0.0);
        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        for i in 1:length(globalNumW1)
            @views @. w11+=wval[globalNumW1[i]]*jphiWn1[1,:,i];
            @views @. w12+=wval[globalNumW1[i]]*jphiWn1[2,:,i];
            @views @. w21+=wval[globalNumW2[i]]*jphiWn2[1,:,i];
            @views @. w22+=wval[globalNumW2[i]]*jphiWn2[2,:,i];
        end

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        fill!(lM11,0.0);
        fill!(lM12,0.0);
        fill!(lM21,0.0);
        fill!(lM22,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk
                    lM11[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ1[r]^2*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j])*(w11[r]*jphiTn1[1,r,i]+w12[r]*jphiTn1[2,r,i]);
                    lM12[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ2[r]^2*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j])*(w21[r]*jphiTn1[1,r,i]+w22[r]*jphiTn1[2,r,i]);
                    lM21[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ1[r]^2*(n[1]*phiFn1[1,r,j]+n[2]*phiFn1[2,r,j])*(w11[r]*jphiTn2[1,r,i]+w12[r]*jphiTn2[2,r,i]);
                    lM22[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ2[r]^2*(n[1]*phiFn2[1,r,j]+n[2]*phiFn2[2,r,j])*(w21[r]*jphiTn2[1,r,i]+w22[r]*jphiTn2[2,r,i]);
                end
            end
        end


        l2g!(globalNumF1,degFF,inc1);
        l2g!(globalNumT1,degFT,inc1);
        l2g!(globalNumF2,degFF,inc2);
        l2g!(globalNumT2,degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];
                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end

        end

    end
    return nothing;
end

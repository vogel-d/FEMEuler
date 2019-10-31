function discGalerkinEdges!(M::Array{Float64,2},
                            degFT::degF{1},phiT::Array{Array{Float64,2},1}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{1},phiW::Array{Array{Float64,2},1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64)



    nT=size(phiT,2);
    nF=size(phiF,2);
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

        phiFn1=@views phiFtrans[eT1];
        phiTn1=@views phiTtrans[eT1];
        phiWn1=@views phiWtrans[eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[eT2];
        phiTn2=@views phiTtrans[eT2];
        phiWn2=@views phiWtrans[eT2];
        kubPn2=@views nquadPoints[eT2];

        fill!(w1,0.0);
        fill!(w2,0.0);
        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        for i in 1:length(globalNumW1)
            @. w1+=wval[globalNumW1[i]]*phiWn1[i];
            @. w2+=wval[globalNumW2[i]]*phiWn2[i];
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
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
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
                            degFT::degF{1},phiT::Array{Array{Float64,2},1}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::Array{Float64,1}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{1},phiW::Array{Array{Float64,2},1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64)



    nT=length(phiT);
    nF=size(phiF,2);
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

        phiFn1=@views phiFtrans[eT1];
        phiTn1=@views phiTtrans[eT1];
        phiWn1=@views phiWtrans[eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[eT2];
        phiTn2=@views phiTtrans[eT2];
        phiWn2=@views phiWtrans[eT2];
        kubPn2=@views nquadPoints[eT2];

        fill!(w1,0.0);
        fill!(w2,0.0);
        l2g!(globalNumW1,degFW,inc1);
        l2g!(globalNumW2,degFW,inc2);
        for i in 1:length(globalNumW1)
            @. w1+=wval[globalNumW1[i]]*phiWn1[i];
            @. w2+=wval[globalNumW2[i]]*phiWn2[i];
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
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
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
                            degFT::degF{2},phiT::Array{Array{Float64,2},2}, phiTtrans::Array{Array{Array{Float64,1},2},1}, globalNumT1::Array{Int64,1}, globalNumT2::Array{Int64,1},
                            degFF::degF{2},phiF::Array{Array{Float64,2},2}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64}, globalNumF1::Array{Int64,1}, globalNumF2::Array{Int64,1},
                            degFW::degF{2},phiW::Array{Array{Float64,2},2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1}, globalNumW1::Array{Int64,1}, globalNumW2::Array{Int64,1},
                            m::mesh, quadWeights::Array{Float64,1}, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},gamma::Float64, coord::Array{Float64,2})



    nT=size(phiT,2);
    nF=size(phiF,2);
    sk=length(quadWeights)

    J1=initPhi((2,2),sk);
    ddJ1=Array{Float64,1}(undef,sk);
    jphiWn1=initPhi(size(phiW),sk)
    jphiTn1=initPhi(size(phiT),sk)

    J2=initPhi((2,2),sk);
    ddJ2=Array{Float64,1}(undef,sk);
    jphiWn2=initPhi(size(phiW),sk)
    jphiTn2=initPhi(size(phiT),sk);

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

        phiFn1=@views phiFtrans[eT1];
        phiTn1=@views phiTtrans[eT1];
        phiWn1=@views phiWtrans[eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[eT2];
        phiTn2=@views phiTtrans[eT2];
        phiWn2=@views phiWtrans[eT2];
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
            @. w11+=wval[globalNumW1[i]]*jphiWn1[1,i];
            @. w12+=wval[globalNumW1[i]]*jphiWn1[2,i];
            @. w21+=wval[globalNumW2[i]]*jphiWn2[1,i];
            @. w22+=wval[globalNumW2[i]]*jphiWn2[2,i];
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
                for r in 1:length(quadWeights)
                    lM11[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ1[r]^2*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r])*(w11[r]*jphiTn1[1,i][r]+w12[r]*jphiTn1[2,i][r]);
                    lM12[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ2[r]^2*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r])*(w21[r]*jphiTn1[1,i][r]+w22[r]*jphiTn1[2,i][r]);
                    lM21[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ1[r]^2*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r])*(w11[r]*jphiTn2[1,i][r]+w12[r]*jphiTn2[2,i][r]);
                    lM22[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ2[r]^2*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r])*(w21[r]*jphiTn2[1,i][r]+w22[r]*jphiTn2[2,i][r]);
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

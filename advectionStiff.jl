function advectionStiff(degFT::degF{1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2}, phiFtrans::Array{Array{Array{Float64,1},2},1},  fval::SparseVector{Float64,Int64},
                        degFW::degF{1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})

    phiF=@views degFF.phi;
    phiW=@views degFW.phi;
    phiT=@views degFT.phi;
    dphiF=@views degFF.divphi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);
    nT=length(phiT);
    nF=size(phiF,2);

    M=zeros(size(degFT.coordinates,2),1);
    for k in 1:m.topology.size[3]

        globalNumW=@views l2g(degFW,k);

        w=zeros(sk);
        for i in 1:length(globalNumW)
            w+=wval[globalNumW[i]]*phiW[i];
        end

        gradw1=zeros(sk);
        gradw2=zeros(sk);
        for r in 1:sk[2]
            for l in 1:sk[1]
                for j in 1:size(gradphiW,2)
                    gradw1[l,r]+=wval[globalNumW[j]]*gradphiW[1,j][l,r]
                    gradw2[l,r]+=wval[globalNumW[j]]*gradphiW[2,j][l,r];
                end
            end
        end

        globalNumF=@views l2g(degFF,k);
        globalNumT=@views l2g(degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        z+=fval[gj]*kubWeights[l,r]*phiT[i][l,r]*(w[l,r]*dphiF[j][l,r]+(gradw1[l,r]*phiF[1,j][l,r]+gradw2[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
            M[gi]-=z;
        end
    end

    quadPoints, quadWeights=getQuad(2*sk[1]-1);

    J1=Array{Array{Float64,1},2}(undef,2,2);
    J2=Array{Array{Float64,1},2}(undef,2,2);
    dJ1=Array{Float64,1}(undef,sk[1]);
    dJ2=Array{Float64,1}(undef,sk[1]);

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

        jacobi!(J1, dJ1,m,inc1,kubPn1);
        jacobi!(J2, dJ2,m,inc2,kubPn2);

        w1=zeros(sk[1])
        w2=zeros(sk[1])
        globalNumW1=@views l2g(degFW,inc1);
        globalNumW2=@views l2g(degFW,inc2);
        for i in 1:length(globalNumW1)
            w1+=wval[globalNumW1[i]]*phiWn1[i];
            w2+=wval[globalNumW2[i]]*phiWn2[i];
        end

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        lM11=zeros(nT,nF);
        lM12=zeros(nT,nF);
        lM21=zeros(nT,nF);
        lM22=zeros(nT,nF);

        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                end
            end
        end
        globalNumF1=@views l2g(degFF,inc1);
        globalNumT1=@views l2g(degFT,inc1);
        globalNumF2=@views l2g(degFF,inc2);
        globalNumT2=@views l2g(degFT,inc2);

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
    return M[1:degFT.num]
end

function advectionStiff(degFT::degF{2}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2}, phiFtrans::Array{Array{Array{Float64,1},2},1},  fval::SparseVector{Float64,Int64},
                        degFW::degF{2}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})


    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    phiW=@views degFW.phi;
    phiT=@views degFT.phi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);
    nT=size(phiT,2);
    nF=size(phiF,2);

    M=zeros(size(degFT.coordinates,2),1);

    ddJ=Array{Float64,2}(undef,sk);
    jphiT=initPhi(size(phiT),sk);

    for k in 1:m.topology.size[3]
        jacobi!(ddJ,jphiT,m,k,kubPoints,phiT);
        globalNumW=@views l2g(degFW,k);

        w1=zeros(sk);
        w2=zeros(sk);
        for i in 1:length(globalNumW)
            w1+=wval[globalNumW[i]]*phiW[1,i];
            w2+=wval[globalNumW[i]]*phiW[2,i];
        end

        gradw11=zeros(sk);
        gradw12=zeros(sk);
        gradw21=zeros(sk);
        gradw22=zeros(sk);
        zg=0;
        for i in 1:size(phiW,2)
            gradw11+=wval[globalNumW[i]]*gradphiW[1,1+zg];
            gradw12+=wval[globalNumW[i]]*gradphiW[1,2+zg];
            gradw21+=wval[globalNumW[i]]*gradphiW[2,1+zg];
            gradw22+=wval[globalNumW[i]]*gradphiW[2,2+zg];
            zg+=2;
        end

        globalNumF=@views l2g(degFF,k);
        globalNumT=@views l2g(degFT,k);
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        z+=fval[gj]*kubWeights[l,r]*ddJ[l,r]^2*(jphiT[1,i][l,r]*(dphiF[j][l,r]*w1[l,r]+gradw11[l,r]*phiF[1,j][l,r]+gradw12[l,r]*phiF[2,j][l,r])+jphiT[2,i][l,r]*(dphiF[j][l,r]*w2[l,r]+gradw21[l,r]*phiF[1,j][l,r]+gradw22[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
            M[gi]-=z;
        end
    end
    quadPoints, quadWeights=getQuad(2*sk[1]-1);

    J1=Array{Array{Float64,1},2}(undef,2,2);
    ddJ1=Array{Float64,1}(undef,sk[1]);
    jphiWn1=initPhi(size(phiW),sk[1])
    jphiTn1=initPhi(size(phiT),sk[1])

    J2=Array{Array{Float64,1},2}(undef,2,2);
    ddJ2=Array{Float64,1}(undef,sk[1]);
    jphiWn2=initPhi(size(phiW),sk[1])
    jphiTn2=initPhi(size(phiT),sk[1]);

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

        jacobi!(J1,ddJ1,jphiWn1,jphiTn1,m,inc1,kubPn1, phiWn1, phiTn1);
        jacobi!(J2,ddJ2,jphiWn2,jphiTn2,m,inc2,kubPn2, phiWn2, phiTn2);


        w11=zeros(sk[1])
        w12=zeros(sk[1])
        w21=zeros(sk[1])
        w22=zeros(sk[1])
        globalNumW1=@views l2g(degFW,inc1);
        globalNumW2=@views l2g(degFW,inc2);
        for i in 1:length(globalNumW1)
            w11+=wval[globalNumW1[i]]*jphiWn1[1,i];
            w12+=wval[globalNumW1[i]]*jphiWn1[2,i];
            w21+=wval[globalNumW2[i]]*jphiWn2[1,i];
            w22+=wval[globalNumW2[i]]*jphiWn2[2,i];
        end

        s=0.0;
        for i in globv
            s+=fval[i];
        end
        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        lM11=zeros(nT,nF);
        lM12=zeros(nT,nF);
        lM21=zeros(nT,nF);
        lM22=zeros(nT,nF);

        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    lM11[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ1[r]^2*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r])*(w11[r]*jphiTn1[1,i][r]+w12[r]*jphiTn1[2,i][r]);
                    lM12[i,j]+=le^2*quadWeights[r]*ddJ1[r]*ddJ2[r]^2*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r])*(w21[r]*jphiTn1[1,i][r]+w22[r]*jphiTn1[2,i][r]);
                    lM21[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ1[r]^2*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r])*(w11[r]*jphiTn2[1,i][r]+w12[r]*jphiTn2[2,i][r]);
                    lM22[i,j]+=le^2*quadWeights[r]*ddJ2[r]*ddJ2[r]^2*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r])*(w21[r]*jphiTn2[1,i][r]+w22[r]*jphiTn2[2,i][r]);
                end
            end
        end


        phiFdeg1=@views l2g(degFF,inc1);
        globalNumT1=@views l2g(degFT,inc1);
        phiFdeg2=@views l2g(degFF,inc2);
        globalNumT2=@views l2g(degFT,inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(phiFdeg2)
                gj1=phiFdeg1[j];
                gj2=phiFdeg2[j];
                M[gi1]+=(+0.5-gammaLoc)*lM11[i,j]*fval[gj1];
                M[gi1]+=(-0.5+gammaLoc)*lM12[i,j]*fval[gj2];
                M[gi2]+=(+0.5+gammaLoc)*lM21[i,j]*fval[gj1];
                M[gi2]+=(-0.5-gammaLoc)*lM22[i,j]*fval[gj2];
            end

        end
    end
    return M[1:degFT.num];
end

function advectionStiff(degFT::degF{1}, phiTtrans::Array{Array{Array{AbstractFloat,1},2},1},
                        degFF::degF{2}, phiFtrans::Array{Array{Array{AbstractFloat,1},2},1},  fval::SparseVector{AbstractFloat,Int},
                        degFW::degF{1}, phiWtrans::Array{Array{Array{AbstractFloat,1},2},1}, wval::Array{AbstractFloat,1},
                        gamma::AbstractFloat,m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2},
                        nquadPoints::Array{Array{AbstractFloat,2},1}, edgeData::Array{Array{Int,1},1})

    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    phiW=@views degFW.phi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);

    globalNumT1=Array{Int,1}(undef,length(phiT));
    globalNumF1=Array{Int,1}(undef,size(phiF,2));
    globalNumW1=Array{Int,1}(undef,length(phiW));

    globalNumT2=Array{Int,1}(undef,length(phiT));
    globalNumF2=Array{Int,1}(undef,size(phiF,2));
    globalNumW2=Array{Int,1}(undef,length(phiW));

    M=zeros(degFT.numB,1);

    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights)

    quadPoints, quadWeights=getQuad(2*sk[1]-1);
    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma)
    return M[1:degFT.num]
end

function advectionStiff(degFT::degF{2}, phiTtrans::Array{Array{Array{AbstractFloat,1},2},1},
                        degFF::degF{2}, phiFtrans::Array{Array{Array{AbstractFloat,1},2},1},  fval::SparseVector{AbstractFloat,Int},
                        degFW::degF{2}, phiWtrans::Array{Array{Array{AbstractFloat,1},2},1}, wval::Array{AbstractFloat,1},
                        gamma::AbstractFloat,m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2},
                        nquadPoints::Array{Array{AbstractFloat,2},1}, edgeData::Array{Array{Int,1},1})


    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    gradphiW=@views degFW.gradphi;
    phiW=@views degFW.phi;

    sk=size(kubWeights);

    quadPoints, quadWeights=getQuad(2*sk[1]-1);
    coord=Array{AbstractFloat,2}(undef,2,m.meshType);

    globalNumT1=Array{Int,1}(undef,size(phiT,2));
    globalNumF1=Array{Int,1}(undef,size(phiF,2));
    globalNumW1=Array{Int,1}(undef,size(phiW,2));

    globalNumT2=Array{Int,1}(undef,size(phiT,2));
    globalNumF2=Array{Int,1}(undef,size(phiF,2));
    globalNumW2=Array{Int,1}(undef,size(phiW,2));

    M=zeros(degFT.numB,1);
    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)

    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma,coord)

    return M[1:degFT.num];
end

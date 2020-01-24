function advectionStiffMatrix(degFT::degF{1}, phiTtrans::Array{Array{Array{AbstractFloat,1},2},1},
                              degFF::degF{2}, phiFtrans::Array{Array{Array{AbstractFloat,1},2},1},fval::Array{AbstractFloat,1},
                              degFW::degF{1}, phiWtrans::Array{Array{Array{AbstractFloat,1},2},1},wval::Array{AbstractFloat,1},
                              gamma::AbstractFloat, m::mesh, kubPoints::Array{AbstractFloat,2}, kubWeights::Array{AbstractFloat,2},
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

    rows=Int[];
    cols=Int[];
    vals=AbstractFloat[];


    discGalerkinCells!(rows, cols, vals,
                       degFT,phiT, globalNumT1,
                       degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights)


    quadPoints, quadWeights=getQuad(2*sk[1]-1);

    discGalerkinEdges!(rows, cols, vals,
                       degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma)


    return sparse(rows,cols,vals)[1:degFT.num,:]
end

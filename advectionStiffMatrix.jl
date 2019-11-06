function advectionStiffMatrix(degFT::degF{3}, phiTtrans::Array{Float64,4},
                              degFF::degF{4}, phiFtrans::Array{Float64,4},fval::Array{Float64,1},
                              degFW::degF{3}, phiWtrans::Array{Float64,4},wval::Array{Float64,1},
                              gamma::Float64, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                              nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})


    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    phiW=@views degFW.phi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);

    globalNumT1=Array{Int64,1}(undef,size(phiT,3));
    globalNumF1=Array{Int64,1}(undef,size(phiF,4));
    globalNumW1=Array{Int64,1}(undef,size(phiW,3));

    globalNumT2=Array{Int64,1}(undef,size(phiT,3));
    globalNumF2=Array{Int64,1}(undef,size(phiF,4));
    globalNumW2=Array{Int64,1}(undef,size(phiW,3));

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];


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

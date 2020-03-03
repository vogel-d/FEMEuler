function advectionStiffMatrix(degFT::degF{1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                              degFF::degF{2}, phiFtrans::Array{Array{Array{Float64,1},2},1},fval::Array{Float64,1},
                              degFW::degF{1}, phiWtrans::Array{Array{Array{Float64,1},2},1},wval::Array{Float64,1},
                              gamma::Float64, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                              nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})


    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    phiW=@views degFW.phi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);

    globalNumT1=Array{Int64,1}(undef,length(phiT));
    globalNumF1=Array{Int64,1}(undef,size(phiF,2));
    globalNumW1=Array{Int64,1}(undef,length(phiW));

    globalNumT2=Array{Int64,1}(undef,length(phiT));
    globalNumF2=Array{Int64,1}(undef,size(phiF,2));
    globalNumW2=Array{Int64,1}(undef,length(phiW));

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    discGalerkinCells!(rows, cols, vals,
                       degFT,phiT, globalNumT1,
                       degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)


    quadPoints, quadWeights=getQuad(2*sk[1]-1);

    discGalerkinEdges!(rows, cols, vals,
                       degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma)


    return sparse(rows,cols,vals)[1:degFT.num,:]
end

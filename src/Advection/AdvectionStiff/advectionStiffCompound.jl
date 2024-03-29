function advectionStiff(degFT::degF{1,:H1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64},
                        degFW::degF{1,:H1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, data::femData)

    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    phiW=@views degFW.phi;
    gradphiW=@views degFW.gradphi;

    sk=size(kubWeights);
    nCompoundPhiT=data.compoundData.nCompoundPhi[degFT.femType];
    nCompoundPhiF=data.compoundData.nCompoundPhi[degFF.femType];
    nCompoundPhiW=data.compoundData.nCompoundPhi[degFW.femType];

    globalNumT1=Array{Int64,1}(undef,nCompoundPhiT);
    globalNumF1=Array{Int64,1}(undef,nCompoundPhiF);
    globalNumW1=Array{Int64,1}(undef,nCompoundPhiW);

    globalNumT2=Array{Int64,1}(undef,nCompoundPhiT);
    globalNumF2=Array{Int64,1}(undef,nCompoundPhiF);
    globalNumW2=Array{Int64,1}(undef,nCompoundPhiW);

    M=zeros(degFT.numB,1);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord, data.compoundData)

    quadWeights=data.compoundData.quadWeights;

    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, data.edgeData, gamma, data.compoundData)

    return M[1:degFT.num]
end

function advectionStiff(degFT::degF{2,:H1div}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1},  fval::SparseVector{Float64,Int64},
                        degFW::degF{2,S} where S, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, data::femData)

    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    gradphiW=@views degFW.gradphi;
    phiW=@views degFW.phi;

    sk=size(kubWeights);

    quadWeights=data.compoundData.quadWeights;

    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    nCompoundPhiT=data.compoundData.nCompoundPhi[degFT.femType];
    nCompoundPhiF=data.compoundData.nCompoundPhi[degFF.femType];
    nCompoundPhiW=data.compoundData.nCompoundPhi[degFW.femType];

    globalNumT1=Array{Int64,1}(undef,nCompoundPhiT);
    globalNumF1=Array{Int64,1}(undef,nCompoundPhiF);
    globalNumW1=Array{Int64,1}(undef,nCompoundPhiW);

    globalNumT2=Array{Int64,1}(undef,nCompoundPhiT);
    globalNumF2=Array{Int64,1}(undef,nCompoundPhiF);
    globalNumW2=Array{Int64,1}(undef,nCompoundPhiW);

    M=zeros(degFT.numB,1);
    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord, data.compoundData)

    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, data.edgeData, gamma, coord, data.compoundData)

    return M[1:degFT.num];
end

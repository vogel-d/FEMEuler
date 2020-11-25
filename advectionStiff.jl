
function advectionStiff(degFT::degF{1,:H1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64},
                        degFW::degF{1,:H1}, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},
                        compoundData::compoundData)

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

    M=zeros(degFT.numB,1);
    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)

    if sk[1]==1
        quadPoints, quadWeights=getQuad(sk[2]+1); #TODO: find better way to handle g (argument of getQuad)
    else
        quadPoints, quadWeights=getQuad(2*sk[1]-1);
    end

    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma, coord)
    return M[1:degFT.num]
end

function advectionStiff(degFT::degF{2,:H1div}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1},  fval::SparseVector{Float64,Int64},
                        degFW::degF{2,S} where S, phiWtrans::Array{Array{Array{Float64,1},2},1}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1},
                        compoundData::compoundData)

    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    gradphiW=@views degFW.gradphi;
    phiW=@views degFW.phi;

    sk=size(kubWeights);

    if sk[1]==1
        quadPoints, quadWeights=getQuad(sk[2]+1); #TODO: find better way to handle g (argument of getQuad)
    else
        quadPoints, quadWeights=getQuad(2*sk[1]-1);
    end


    coord=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT1=Array{Int64,1}(undef,size(phiT,2));
    globalNumF1=Array{Int64,1}(undef,size(phiF,2));
    globalNumW1=Array{Int64,1}(undef,size(phiW,2));

    globalNumT2=Array{Int64,1}(undef,size(phiT,2));
    globalNumF2=Array{Int64,1}(undef,size(phiF,2));
    globalNumW2=Array{Int64,1}(undef,size(phiW,2));

    M=zeros(degFT.numB,1);
    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)

    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData, gamma, coord)

    return M[1:degFT.num];
end

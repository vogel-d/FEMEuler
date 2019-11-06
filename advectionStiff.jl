function advectionStiff(degFT::degF{3}, phiTtrans::Array{Float64,4},
                        degFF::degF{4}, phiFtrans::Array{Float64,4},  fval::SparseVector{Float64,Int64},
                        degFW::degF{3}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
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
    globalNumW2=Array{Int64,1}(undef,size(phiT,3));

    M=zeros(size(degFT.coordinates,2),1);

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

function advectionStiff(degFT::degF{4}, phiTtrans::Array{Float64,4},
                        degFF::degF{4}, phiFtrans::Array{Float64,4},  fval::SparseVector{Float64,Int64},
                        degFW::degF{4}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})


    phiT=@views degFT.phi;
    phiF=@views degFF.phi;
    dphiF=@views degFF.divphi;
    gradphiW=@views degFW.gradphi;
    phiW=@views degFW.phi;

    sk=size(kubWeights);

    quadPoints, quadWeights=getQuad(2*sk[1]-1);
    coord=Array{Float64,2}(undef,2,m.meshType);

    globalNumT1=Array{Int64,1}(undef,size(phiT,4));
    globalNumF1=Array{Int64,1}(undef,size(phiF,4));
    globalNumW1=Array{Int64,1}(undef,size(phiW,4));

    globalNumT2=Array{Int64,1}(undef,size(phiT,4));
    globalNumF2=Array{Int64,1}(undef,size(phiF,4));
    globalNumW2=Array{Int64,1}(undef,size(phiW,4));

    M=zeros(size(degFT.coordinates,2),1);
    #println("cells:")
    discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)
    #println("edges:")
    discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma,coord)

    return M[1:degFT.num];
end

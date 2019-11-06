function advectionStiff(degFT::degF{2}, phiTtrans::Array{Float64,4},
                        degFF::degF{2}, phiFtrans::Array{Float64,4},  fval::SparseVector{Float64,Int64},
                        degFW::degF{2}, phiWtrans::Array{Float64,4}, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})


    phiT=hdimarray(degFT.phi);
    phiF=hdimarray(degFF.phi);
    dphiF=hdimarray(degFF.divphi);
    gradphiW=hdimarray(degFW.gradphi);
    phiW=hdimarray(degFW.phi);

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
    println("3D:")
    @time discGalerkinCells!(M,degFT,phiT, globalNumT1, degFF,phiF, dphiF, fval, globalNumF1,
                       degFW, phiW, gradphiW, wval, globalNumW1,
                       m, kubPoints, kubWeights, coord)


    @time discGalerkinEdges!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       degFW,phiW, phiWtrans, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma,coord)

    return M[1:degFT.num];
end

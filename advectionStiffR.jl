function advectionStiffR(degFT::degF{1,:H1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1}, fval::SparseVector{Float64,Int64},
                        recoverySpace::Symbol, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})

    phiT=@views degFT.phi;
    gradphiT=@views degFT.gradphi;
    phiF=@views degFF.phi;
    phiW=getPhiRecovery(Val(recoverySpace));

    sk=size(kubWeights);

    globalNumT1=Array{Int64,1}(undef,length(phiT));
    globalNumF1=Array{Int64,1}(undef,size(phiF,2));
    globalNumW1=UnitRange{Int64}(1,length(phiW));

    globalNumT2=Array{Int64,1}(undef,length(phiT));
    globalNumF2=Array{Int64,1}(undef,size(phiF,2));
    globalNumW2=UnitRange{Int64}(1,length(phiW));

    M=zeros(degFT.numB,1);
    coord1=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    coord2=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    discGalerkinCellsR!(M,degFT,gradphiT,globalNumT1,
                       degFF,phiF,fval,globalNumF1,
                       phiW, wval,
                       m, kubPoints, kubWeights, coord1)
    if sk[1]==1
        quadPoints, quadWeights=getQuad(sk[2]+1); #TODO: find better way to handle g (argument of getQuad)
    else
        quadPoints, quadWeights=getQuad(2*sk[1]-1);
    end

    discGalerkinEdgesR!(M,degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       phiW, wval, globalNumW1, globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma, coord1, coord2)
    return M[1:degFT.num]
end

function advectionStiffMatrixR(degFT::degF{1,:H1}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                              degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1},fval::Array{Float64,1},
                              recoverySpace::Symbol, wval::Array{Float64,1},
                              gamma::Float64, m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                              nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})

    phiT=@views degFT.phi;
    gradphiT=@views degFT.gradphi;
    phiF=@views degFF.phi;
    phiW=getPhiRecovery(Val(recoverySpace));

    sk=size(kubWeights);

    globalNumT1=Array{Int64,1}(undef,length(phiT));
    globalNumF1=Array{Int64,1}(undef,size(phiF,2));
    globalNumW1=UnitRange{Int64}(1,length(phiW));

    globalNumT2=Array{Int64,1}(undef,length(phiT));
    globalNumF2=Array{Int64,1}(undef,size(phiF,2));
    globalNumW2=UnitRange{Int64}(1,length(phiW));

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    coord1=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    coord2=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    discGalerkinCellsR!(rows, cols, vals,
                       degFT, gradphiT, globalNumT1,
                       degFF, phiF, fval, globalNumF1,
                       phiW, wval,
                       m, kubPoints, kubWeights, coord1)

    if sk[1]==1
        quadPoints, quadWeights=getQuad(sk[2]+1); #TODO: find better way to handle g (argument of getQuad)
    else
        quadPoints, quadWeights=getQuad(2*sk[1]-1);
    end

    discGalerkinEdgesR!(rows, cols, vals,
                       degFT,phiT, phiTtrans,globalNumT1, globalNumT2,
                       degFF,phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       phiW, wval, globalNumW1,globalNumW2,
                       m, quadWeights, nquadPoints, edgeData,gamma, coord1, coord2)

    return sparse(rows,cols,vals)[1:degFT.num,:]
end

function advectionStiffR(degFT::degF{2,:H1div}, phiTtrans::Array{Array{Array{Float64,1},2},1},
                        degFF::degF{2,:H1div}, phiFtrans::Array{Array{Array{Float64,1},2},1},  fval::SparseVector{Float64,Int64},
                        recoverySpace::Symbol, wval::Array{Float64,1},
                        gamma::Float64,m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2},
                        nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})

    phiT=@views degFT.phi;
    gradphiT=@views degFT.gradphi;
    phiF=@views degFF.phi;
    phiW=getPhiRecovery(Val(recoverySpace));

    sk=size(kubWeights);

    if sk[1]==1
        quadPoints, quadWeights=getQuad(sk[2]+1); #TODO: find better way to handle g (argument of getQuad)
    else
        quadPoints, quadWeights=getQuad(2*sk[1]-1);
    end


    coord1=Array{Float64,2}(undef,m.geometry.dim,m.meshType);
    coord2=Array{Float64,2}(undef,m.geometry.dim,m.meshType);

    globalNumT1=Array{Int64,1}(undef,size(phiT,2));
    globalNumF1=Array{Int64,1}(undef,size(phiF,2));
    globalNumW1=UnitRange{Int64}(1,size(phiW,2));

    globalNumT2=Array{Int64,1}(undef,size(phiT,2));
    globalNumF2=Array{Int64,1}(undef,size(phiF,2));
    globalNumW2=UnitRange{Int64}(1,size(phiW,2));

    M=zeros(degFT.numB,1);
    discGalerkinCellsR!(M,degFT,gradphiT,globalNumT1,
                       degFF,phiF,fval,globalNumF1,
                       phiW,wval,globalNumW1,
                       m, kubPoints, kubWeights, coord1)

    discGalerkinEdgesR!(M, degFT, phiT, phiTtrans, globalNumT1, globalNumT2,
                       degFF, phiF, phiFtrans, fval, globalNumF1, globalNumF2,
                       phiW, wval, globalNumW1, globalNumW2,
                       m, quadWeights, nquadPoints, edgeData, gamma, coord1, coord2)

    return M[1:degFT.num];
end

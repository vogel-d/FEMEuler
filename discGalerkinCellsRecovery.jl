function discGalerkinCellsR!(M::Array{Float64,2},
                            degFT::degF{1,:H1}, gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div},phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            phiW::Array{Function,1}, wval::Array{Float64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    nW=length(phiW)

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    zf=1
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        n=transformation(m,coord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        fill!(w,0.0);
        for i in 1:nW
            for r in 1:size(kubWeights,2)
                for l in 1:size(kubWeights,1)
                    w[l,r]+=wval[zf]*phiW[i](transformRecoveryCoord(n,t1,t2,transformation(m,coord,kubPoints[1,l],kubPoints[2,r])));
                end
            end
            zf+=1
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        phiFgradphiT=0.0;
                        for d in 1:m.topology.dim
                            phiFgradphiT+=phiF[d,j][l,r]*gradphiT[d,i][l,r]
                        end
                        z+=fval[gj]*kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*w[l,r]*phiFgradphiT;
                    end
                end
            end
            M[gi]+=z;
        end
    end

    return nothing;
end

function discGalerkinCellsR!(rows::Array{Int64,1}, cols::Array{Int64,1}, vals::Array{Float64,1},
                            degFT::degF{1,:H1}, gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div}, phiF::Array{Array{Float64,2},2}, fval::Array{Float64,1}, globalNumF::Array{Int64,1},
                            phiW::Array{Function,1}, wval::Array{Float64,1},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    nW=length(phiW)
    nT=length(globalNumT);
    nF=length(globalNumF);

    w=zeros(sk);

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    lM=zeros(nT,nF);
    zf=1
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        n=transformation(m,coord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        fill!(w,0.0);
        for i in 1:nW
            for r in 1:size(kubWeights,2)
                for l in 1:size(kubWeights,1)
                    w[l,r]+=wval[zf]*phiW[i](transformRecoveryCoord(n,t1,t2,transformation(m,coord,kubPoints[1,l],kubPoints[2,r])));
                end
            end
            zf+=1
        end


        fill!(lM,0.0);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        phiFgradphiT=0.0;
                        for d in 1:m.topology.dim
                            phiFgradphiT+=phiF[d,j][l,r]*gradphiT[d,i][l,r]
                        end
                        lM[i,j]+=kubWeights[l,r]*(abs(dJ[l,r])/dJ[l,r])*w[l,r]*phiFgradphiT;
                    end
                end
            end
        end
        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            for j in 1:length(globalNumF)
                push!(rows,gi);
                push!(cols,globalNumF[j]);
                push!(vals,lM[i,j]);
            end
        end
    end

    return nothing;
end


function discGalerkinCellsR!(M::Array{Float64,2},
                            degFT::degF{2,:H1div}, gradphiT::Array{Array{Float64,2},2}, globalNumT::Array{Int64,1},
                            degFF::degF{2,:H1div}, phiF::Array{Array{Float64,2},2}, fval::SparseVector{Float64,Int64}, globalNumF::Array{Int64,1},
                            phiW::Array{Function,2}, wval::Array{Float64,1}, globalNumW::UnitRange{Int64},
                            m::mesh, kubPoints::Array{Float64,2}, kubWeights::Array{Float64,2}, coord::Array{Float64,2})

    sk=size(kubWeights);
    nW=size(phiW,2)

    J=initJacobi((m.geometry.dim,m.topology.dim),sk);
    dJ=Array{Float64,2}(undef,sk);

    w=[zeros(sk) for d in 1:m.geometry.dim]

    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints,coord);

        n=transformation(m,coord,0.5,0.5);
        t1,t2=getTangentialPlane(n)

        globalNumW=(nW*(k-1)+1):(nW*k)
        for d in 1:m.geometry.dim
            fill!(w[d], 0.0)
            for i in 1:nW
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        w[d][l,r]+=wval[globalNumW[i]]*phiW[d,i](transformRecoveryCoord(n,t1,t2,transformation(m,coord,kubPoints[1,l],kubPoints[2,r])));
                    end
                end
            end
        end

        l2g!(globalNumF,degFF,k);
        l2g!(globalNumT,degFT,k);
        zg=0;
        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            z=0.0;
            for j in 1:length(globalNumF)
                gj=globalNumF[j];
                for r in 1:size(kubWeights,2)
                    for l in 1:size(kubWeights,1)
                        val=0.0;
                        for d in 1:m.geometry.dim
                            jgradphiTphiF=0.0;
                            for dt in 1:m.topology.dim
                                gradphiTphiF=0.0;
                                for dti in 1:m.topology.dim
                                    gradphiTphiF+=gradphiT[dt,zg+dti][l,r]*phiF[dti,j][l,r]
                                end
                                jgradphiTphiF+=J[d,dt][l,r]*gradphiTphiF
                            end
                            val+=jgradphiTphiF*w[d][l,r]
                        end
                        z+=fval[gj]*kubWeights[l,r]*abs(dJ[l,r])/(dJ[l,r]^2)*val
                    end
                end
            end
            M[gi]+=z;
            zg+=2;
        end

    end
    return nothing;
end

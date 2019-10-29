function advectionStiffMatrix(p::femProblem, compT::Symbol, phiTtrans::Array{Array{Array{Float64,1},2},1},
                              compF::Symbol, phiFtrans::Array{Array{Array{Float64,1},2},1},fval::Array{Float64,1},
                              compW::Symbol, phiWtrans::Array{Array{Array{Float64,1},2},1},wval::Array{Float64,1},
                              gamma::Float64, nquadPoints::Array{Array{Float64,2},1}, edgeData::Array{Array{Int64,1},1})

    degF=p.degFBoundary;
    m=p.mesh;

    phiF=@views degF[compF].phi;
    phiW=@views degF[compW].phi;
    phiT=@views degF[compT].phi;
    dphiF=@views degF[compF].divphi;
    gradphiW=@views degF[compW].gradphi;

    kubPoints=@views p.kubPoints;
    kubWeights=@views p.kubWeights;

    sk=size(kubWeights);
    nT=length(phiT);
    nF=size(phiF,2);

    rows=Int64[];
    cols=Int64[];
    vals=Float64[];

    for k in 1:m.topology.size[3]
        globalNumW=@views l2g(degF[compW],k);

        w=zeros(sk);
        for i in 1:length(globalNumW)
            w+=wval[globalNumW[i]]*phiW[i];
        end

        gradw1=zeros(sk);
        gradw2=zeros(sk);
        for r in 1:sk[2]
            for l in 1:sk[1]
                for j in 1:size(gradphiW,2)
                    gradw1[l,r]+=wval[globalNumW[j]]*gradphiW[1,j][l,r]
                    gradw2[l,r]+=wval[globalNumW[j]]*gradphiW[2,j][l,r];
                end
            end
        end

        lM=zeros(nT,nF);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    for l in 1:sk[1]
                        lM[i,j]+=kubWeights[l,r]*phiT[i][l,r]*(w[l,r]*dphiF[j][l,r]+(gradw1[l,r]*phiF[1,j][l,r]+gradw2[l,r]*phiF[2,j][l,r]));
                    end
                end
            end
        end

        globalNumF=@views l2g(degF[compF],k);
        globalNumT=@views l2g(degF[compT],k);

        for i in 1:length(globalNumT)
            gi=globalNumT[i];
            for j in 1:length(globalNumF)
                push!(rows,gi);
                push!(cols,globalNumF[j]);
                push!(vals,-lM[i,j]);
            end
        end
    end

    quadPoints, quadWeights=getQuad(2*sk[1]-1);

    J1=Array{Array{Float64,1},2}(undef,2,2);
    J2=Array{Array{Float64,1},2}(undef,2,2);
    dJ1=Array{Float64,1}(undef,sk[1]);
    dJ2=Array{Float64,1}(undef,sk[1]);

    z=1;
    for e in 1:length(edgeData[1])
        inc1=edgeData[2][z];
        inc2=edgeData[2][z+1];
        eT1=edgeData[3][z];
        eT2=edgeData[3][z+1];
        z+=2;
        n=@views m.normals[:,eT1];;
        globv=@views edgeData[4][edgeData[5][e]:edgeData[5][e+1]-1];

        phiFn1=@views phiFtrans[eT1];
        phiTn1=@views phiTtrans[eT1];
        wn1=@views phiWtrans[eT1];
        kubPn1=@views nquadPoints[eT1];
        phiFn2=@views phiFtrans[eT2];
        phiTn2=@views phiTtrans[eT2];
        phiWn2=@views phiWtrans[eT2];
        kubPn2=@views nquadPoints[eT2];

        jacobi!(J1, dJ1,m,inc1,kubPn1);
        jacobi!(J2, dJ2,m,inc2,kubPn2);

        w1=zeros(sk[1])
        w2=zeros(sk[1])
        globalNumW1=@views l2g(degF[compW],inc1);
        globalNumW2=@views l2g(degF[compW],inc2);
        for i in 1:length(globalNumW1)
            w1+=wval[globalNumW1[i]]*wn1[i];
            w2+=wval[globalNumW2[i]]*phiWn2[i];
        end

        s=0.0;
        for i in globv
            i>length(fval) && continue;
            s+=fval[i];
        end

        s<=0.0 ? gammaLoc=-gamma : gammaLoc=gamma;

        lM11=zeros(nT,nF);
        lM12=zeros(nT,nF);
        lM21=zeros(nT,nF);
        lM22=zeros(nT,nF);
        for j in 1:nF
            for i in 1:nT
                for r in 1:sk[2]
                    lM11[i,j]+=quadWeights[r]*w1[r]*phiTn1[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM12[i,j]+=quadWeights[r]*w2[r]*phiTn1[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                    lM21[i,j]+=quadWeights[r]*w1[r]*phiTn2[i][r]*(n[1]*phiFn1[1,j][r]+n[2]*phiFn1[2,j][r]);
                    lM22[i,j]+=quadWeights[r]*w2[r]*phiTn2[i][r]*(n[1]*phiFn2[1,j][r]+n[2]*phiFn2[2,j][r]);
                end
            end
        end

        globalNumF1=@views l2g(degF[compF],inc1);
        globalNumT1=@views l2g(degF[compT],inc1);
        globalNumF2=@views l2g(degF[compF],inc2);
        globalNumT2=@views l2g(degF[compT],inc2);

        for i in 1:length(globalNumT1)
            gi1=globalNumT1[i];
            gi2=globalNumT2[i];
            for j in 1:length(globalNumF2)
                gj1=globalNumF1[j];
                gj2=globalNumF2[j];

                push!(rows,gi1);
                push!(cols,gj1);
                push!(vals,(+0.5-gammaLoc)*lM11[i,j]);

                push!(rows,gi1);
                push!(cols,gj2);
                push!(vals,(-0.5+gammaLoc)*lM12[i,j]);

                push!(rows,gi2);
                push!(cols,gj1);
                push!(vals,(+0.5+gammaLoc)*lM21[i,j]);

                push!(rows,gi2);
                push!(cols,gj2);
                push!(vals,(-0.5-gammaLoc)*lM22[i,j]);
            end
        end
    end

    return sparse(rows,cols,vals)[1:degF[compT].num,:]
end

function projectPressure(p::femProblem,compP::Symbol,compRT::Symbol,valRT::Array{Float64,1})
    degF=p.degFBoundary;
    m=p.mesh;
    phiRT=@views degF[compRT].phi;
    phiP=@views degF[compP].phi;
    kubPoints=@views p.kubPoints;
    kubWeights=@views p.kubWeights;
    sk=size(kubWeights);

    J=Array{Array{Float64,2},2}(undef,2,2);
    dJ=Array{Float64,2}(undef,sk);

    gbh=zeros(size(degF[compP].coordinates,2))
    for k in 1:m.topology.size[3]
        jacobi!(J,dJ,m,k,kubPoints);

        globalNumRT=@views l2g(degF[compRT],k);
        globalNumP=@views l2g(degF[compP],k);

        cl=zeros(sk);
        for i in 1:length(globalNumRT)
            cl+=valRT[globalNumRT[i]]*phiRT[1,i];
        end

        Cpd=1004.0;
        Rd=Cpd-717.0; kappa=Rd/Cpd;
        p0=100000.0;
        Pres1=p0*(Rd/p0)^(1.0/(1.0-kappa));
        Pres2=(1.0/(1.0-kappa));
        cl=Pres1*cl.^Pres2;

        for j in 1:length(phiP)
            for r in 1:sk[2]
                for l in 1:sk[1]
                    gbh[globalNumP[j]]+=kubWeights[l,r]*cl[l,r]*phiP[1,j][l,r]*dJ[l,r];
                end
            end
        end
    end
    return p.massMBoundary[compP]\gbh;
end

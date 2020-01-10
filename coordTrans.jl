function coordTrans(mt::Int64, normals::Array{Float64,2}, type::Array{Symbol,1}, n::Int64)
    g=2*n-1;
    quadPoints, quadWeights=getQuad(g);
    sk=length(quadPoints);

    nquadPhi=Dict{Symbol, Array{Array{Array{Float64,1},2},1}}();
    nquadPoints=Array{Array{Float64,2},1}(undef, size(normals,2));

    for i in 1:size(normals,2)
        n=normals[:,i];
        if n==[0.0,-1.0]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[1,:]=quadPoints;
        elseif n==[-1.0,0.0]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[2,:]=quadPoints;
        elseif n==[0.0, 1.0]
            newQuadPoints=ones(2,sk);
            newQuadPoints[1,:]=quadPoints;
        elseif n==[1.0,0.0]
            newQuadPoints=ones(2,sk);
            newQuadPoints[2,:]=quadPoints;
        elseif n==[0.7071067811865475244, 0.7071067811865475244]
            newQuadPoints=zeros(2,sk);
            newQuadPoints[2,:]=quadPoints;
            newQuadPoints[1,:]=1 .- quadPoints;
        end
        nquadPoints[i]=newQuadPoints;
    end

    for k in type
        nquadPhi[k]=Array{Array{Array{Float64,1},2},1}(undef, size(normals,2));
        phi, psize =getElementProperties(k,mt);
        for m in 1:size(normals,2)
            quadPhi=Array{Array{Float64,1},2}(undef,psize[1], psize[2]);
            for n in 1:length(phi)
                quadVal=Array{Float64,1}(undef,sk);
                for i=1:sk
                    quadVal[i]=phi[n](nquadPoints[m][1,i], nquadPoints[m][2,i]);
                end
                quadPhi[n]=quadVal;
            end
            nquadPhi[k][m]=quadPhi;
        end
    end
    return nquadPhi, nquadPoints;
end

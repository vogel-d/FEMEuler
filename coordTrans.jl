function coordTrans(mt::Int64, normals::Array{Float64,2}, type::Array{Symbol,1}, n::Int64)
    quadPoints, quadWeights=getQuad(2*n-1);
    sk=length(quadPoints);

    nquadPhi=Dict{Symbol, Array{Float64,4}}();
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
        phi, psize =getPhi(k);
        nquadPhi[k]=Array{Float64,4}(undef,psize[1],sk,psize[2],size(normals,2));
        for m in 1:size(normals,2)
            for l in 1:psize[2]
                for n in 1:psize[1]
                    for i=1:sk
                        nquadPhi[k][n,i,l,m]=phi[l,n](nquadPoints[m][1,i], nquadPoints[m][2,i]);
                    end
                end
            end
        end
    end
    return nquadPhi, nquadPoints;
end

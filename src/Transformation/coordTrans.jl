function coordTrans(mesh::mesh, normals::Array{Float64,2}, type::Array{Symbol,1}, n::Int64)
    g=2*n-1;
    quadPoints, quadWeights=getQuad(g);
    sk=length(quadPoints);
    mt=mesh.meshType

    nquadPhi=Dict{Symbol, Array{Array{Array{Float64,1},2},1}}();
    nquadPoints=Array{Array{Float64,2},1}(undef, size(normals,2));
    if mt==4
        #Kante 1
        nquadPoints[1]=zeros(2,sk);
        nquadPoints[1][1,:]=quadPoints;
        #nquadPoints[1][1,:]=ones(sk).-quadPoints;

        #Kante 2
        nquadPoints[2]=ones(2,sk);
        nquadPoints[2][2,:]=quadPoints;
        #nquadPoints[2][2,:]=ones(sk).-quadPoints;

        #Kante 3
        nquadPoints[3]=ones(2,sk);
        nquadPoints[3][1,:]=quadPoints;
        #nquadPoints[3][1,:]=ones(sk).-quadPoints;

        #Kante 4
        nquadPoints[4]=zeros(2,sk);
        nquadPoints[4][2,:]=quadPoints;
        #nquadPoints[4][2,:]=ones(sk).-quadPoints;

    elseif mt==3
        #Kante 1
        nquadPoints[1]=zeros(2,sk);
        nquadPoints[1][1,:]=quadPoints;

        #Kante 2
        nquadPoints[2]=zeros(2,sk);
        nquadPoints[2][1,:]=1 .- quadPoints;
        nquadPoints[2][2,:]=quadPoints;

        #Kante 3
        nquadPoints[3]=zeros(2,sk);
        nquadPoints[3][2,:]=quadPoints;

    else
        error("Unbekannter meshType")
    end
    nf=mesh.topology.size[3];
    mcoord=mesh.geometry.coordinates;
    inc=mesh.topology.incidence["20"];
    off=mesh.topology.offset["20"];


    for k in type
        nquadPhi[k]=Array{Array{Array{Float64,1},2},1}(undef, size(normals,2));
        phi, psize =getElementProperties(k,mt,mesh.geometry.dim);
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

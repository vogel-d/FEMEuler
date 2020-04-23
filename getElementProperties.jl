include("getQuadElementProperties.jl")
include("getTriElementProperties.jl")

function getElementProperties(type::Symbol, kubPoints::Array{Float64,2}, mt::Int)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);

        sk=size(kubPoints,2);
        kubPhi=Array{Array{Float64,2},ndims(phi)}(undef,size(phi));
        for k=1:length(phi)
            kubVal=Array{Float64,2}(undef,sk,sk);
            for i=1:sk, j=1:sk
                kubVal[i,j]=phi[k](kubPoints[1,i], kubPoints[2,j]);
                #für Dreiecke werden hier wesentlich mehr (mögl. redundante) Werte als nötig ausgerechnet,
                #nur die Hauptdiagonale ist dann relevant,
                #wird über kubWeights geregelt
            end
            kubPhi[k]=kubVal;
        end

        kubDiv=Array{Array{Float64,2},1}(undef,length(divphi));
        for k in 1:length(divphi)
            kubVal=Array{Float64,2}(undef,sk,sk);
            for i=1:sk, j=1:sk
                kubVal[i,j]=divphi[k](kubPoints[1,i], kubPoints[2,j]);
            end
            kubDiv[k]=kubVal;
        end

        kubGrad=Array{Array{Float64,2},2}(undef,size(gradphi));
        for ki=1:size(gradphi,1), kj=1:size(gradphi,2)
            kubVal=Array{Float64,2}(undef,sk,sk);
            for i=1:sk, j=1:sk
                kubVal[i,j]=gradphi[ki,kj](kubPoints[1,i], kubPoints[2,j]);
            end
            kubGrad[ki,kj]=kubVal;
        end

    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);

        sk=size(kubPoints,2);
        kubPhi=Array{Array{Float64,2},ndims(phi)}(undef,size(phi));
        for k=1:length(phi)
            kubVal=Array{Float64,2}(undef,1,sk);
            for i=1:sk
                kubVal[1,i]=phi[k](kubPoints[1,i], kubPoints[2,i]);
                #für Dreiecke werden hier wesentlich mehr (mögl. redundante) Werte als nötig ausgerechnet,
                #nur die Hauptdiagonale ist dann relevant,
                #wird über kubWeights geregelt
            end
            kubPhi[k]=kubVal;
        end

        kubDiv=Array{Array{Float64,2},1}(undef,length(divphi));
        for k in 1:length(divphi)
            kubVal=Array{Float64,2}(undef,1,sk);
            for i=1:sk
                kubVal[1,i]=divphi[k](kubPoints[1,i], kubPoints[2,i]);
            end
            kubDiv[k]=kubVal;
        end

        kubGrad=Array{Array{Float64,2},2}(undef,size(gradphi));
        for ki=1:size(gradphi,1), kj=1:size(gradphi,2)
            kubVal=Array{Float64,2}(undef,1,sk);
            for i=1:sk
                kubVal[1,i]=gradphi[ki,kj](kubPoints[1,i], kubPoints[2,i]);
            end
            kubGrad[ki,kj]=kubVal;
        end

    end

    return kubPhi, kubDiv,  kubGrad, nFace, nEdge, nVert
end

function getElementProperties(type::Symbol, kubPoints::Array{Float64,2}, m::mesh)

    sk=size(kubPoints,2);
    nf=m.topology.size[3]

    phi, divphi, gradphi, cm, nFace, nEdge, nVert=getRecoveryElementProperties(type);
    if ndims(phi)==1
        nPhi=length(phi)
        kubPhi=Array{Array{Float64,2},1}(undef,nf*nPhi);
    else
        nPhi=size(phi,2)
        kubPhi=Array{Array{Float64,2},2}(undef,size(phi,1),nf*nPhi);
    end
    kubDiv=Array{Array{Float64,2},1}(undef,0);
    nGradPhi=size(gradphi,2)
    kubGrad=Array{Array{Float64,2},2}(undef,3,nf*nGradPhi);

    mcoord=m.geometry.coordinates;
    inc=m.topology.incidence["20"];
    off=m.topology.offset["20"];
    z=0; zg=0;
    for k in 1:nf
        coord=@views mcoord[:,inc[off[k]:off[k+1]-1]]
        mp=transformation(m, coord, 0.5, 0.5);


        for k=1:length(phi)
            kubVal=Array{Float64,2}(undef,sk,sk);
            for i=1:sk, j=1:sk
                kubVal[i,j]=phi[k](transformation(m, coord, kubPoints[1,i], kubPoints[2,j]), mp);
            end
            kubPhi[z+k]=kubVal;
        end


        for ki=1:size(gradphi,1), kj=1:size(gradphi,2)
            kubVal=Array{Float64,2}(undef,sk,sk);
            for i=1:sk, j=1:sk
                kubVal[i,j]=gradphi[ki,kj](transformation(m, coord, kubPoints[1,i], kubPoints[2,j]), mp);
            end
            kubGrad[ki,zg+kj]=kubVal;
        end

        z+=nPhi; zg+=nGradPhi;
    end

    return kubPhi, kubDiv, kubGrad, nFace, nEdge, nVert
end

function getElementProperties(type::Symbol, mt::Int, recovery::Bool=false)
    if recovery
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getRecoveryElementProperties(type);
    elseif mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    ndims(phi)==1 ? s=size(phi') : s=size(phi);
    return phi, s
end

function getElementProperties(mt::Int, type::Symbol, recovery::Bool=false)
    if recovery
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getRecoveryElementProperties(type);
    elseif mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    return cm
end

function getElementProperties(type::Symbol, mt::Int, x, y)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    valPhi=similar(phi,Float64);
    for k=1:length(phi)
        valPhi[k]=phi[k](x,y);
    end
    return valPhi
end

function getElementProperties(type::Symbol, mt::Int, mp::Array{Float64,1}, xyz)
    phi, divphi, gradphi, cm, nFace, nEdge, nVert=getRecoveryElementProperties(type);

    valPhi=similar(phi,Float64);
    for k=1:length(phi)
        valPhi[k]=phi[k](xyz,mp);
    end
    return valPhi
end

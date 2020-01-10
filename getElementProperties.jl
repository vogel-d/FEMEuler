include("getQuadElementProperties.jl")
include("getTriElementProperties.jl")

function getElementProperties(type::Symbol, kubPoints::Array{Float64,2}, mt::Int)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert, comp=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    sk=size(kubPoints,2);
    kubPhi=Array{Array{Float64,2},ndims(phi)}(undef,size(phi));
    for k=1:length(phi)
        kubVal=Array{Float64,2}(undef,sk,sk);
        for i=1:sk, j=1:sk
            kubVal[i,j]=phi[k](kubPoints[1,i], kubPoints[2,j]);
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

    return kubPhi, kubDiv,  kubGrad, comp, nFace, nEdge, nVert
end

function getElementProperties(type::Symbol, mt::Int)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert, comp=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    ndims(phi)==1 ? s=size(phi') : s=size(phi);
    return phi, s
end

function getElementProperties(mt::Int, type::Symbol)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert, comp=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    return cm
end

function getElementProperties(type::Symbol, mt::Int, x, y)
    if mt==4
        phi, divphi, gradphi, cm, nFace, nEdge, nVert, comp=getQuadElementProperties(type);
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    valPhi=similar(phi,Float64);
    for k=1:length(phi)
        valPhi[k]=phi[k](x,y);
    end
    return valPhi
end

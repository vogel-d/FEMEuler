include("getQuadElementProperties.jl")
include("getSphereElementProperties.jl")
include("getTriElementProperties.jl")

function getElementProperties(type::Symbol, kubPoints::Array{Float64,2}, mt::Int, dim::Int)
    if mt==4
        if dim==2
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
        elseif dim==3
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getSphereElementProperties(type);
        end

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

function getElementProperties(type::Symbol, mt::Int, dim::Int)
    if mt==4
        if dim==2
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
        elseif dim==3
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getSphereElementProperties(type);
        end
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    ndims(phi)==1 ? s=size(phi') : s=size(phi);
    return phi, s
end

function getElementProperties(mt::Int, type::Symbol, dim::Int)
    if mt==4
        if dim==2
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
        elseif dim==3
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getSphereElementProperties(type);
        end
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    return cm
end

function getElementProperties(type::Symbol, mt::Int, dim::Int, x, y)
    if mt==4
        if dim==2
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getQuadElementProperties(type);
        elseif dim==3
            phi, divphi, gradphi, cm, nFace, nEdge, nVert=getSphereElementProperties(type);
        end
    elseif mt==3
        phi, divphi, gradphi, cm, nFace, nEdge, nVert=getTriElementProperties(type);
    end
    valPhi=similar(phi,Float64);
    for k=1:length(phi)
        valPhi[k]=phi[k](x,y);
    end
    return valPhi
end

function getElementProperties(type::Symbol, n::Array{Float64,1}, xyz)
    phi=getPhiRecovery(Val(type));
    valPhi=similar(phi,Float64);
    t1,t2=getTangentialPlane(n);
    txyz=transformRecoveryCoord(n,t1,t2,xyz)
    for k=1:length(phi)
        valPhi[k]=phi[k](txyz);
    end
    return valPhi
end
#=
function getElementProperties(type::Symbol, n::Array{Float64,1}, xyz)
    phi=getPhiRecovery([0.0,0.0],Val(type));
    valPhi=similar(phi,Float64);
    n2=copy(n)
    t1,t2=getTangentialPlane(n2);
    txyz=transformRecoveryCoord(n2,t1,t2,xyz)
    for k=1:length(phi)
        valPhi[k]=phi[k](txyz);
    end
    return valPhi
end
=#

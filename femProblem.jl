mutable struct femProblem
    mesh::mesh;
    boundaryValues::Dict{Tuple{Symbol,Symbol}, Array{Float64,1}};
    degFBoundary::Dict{Symbol, degF};
    femType::Dict{Symbol, Array{Symbol,1}};
    edgeData::Array{Array{Int64,1},1};
    compoundData::compoundData;
    solution::Dict{Float64, solution};
    massM::Dict{Symbol, SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}};
    massMBoundary::Dict{Symbol, SuiteSparse.UMFPACK.UmfpackLU{Float64,Int64}};
    stiffM::Dict{Symbol, SparseMatrixCSC{Float64,Int64}};
    recoveryM::Dict{Tuple{Symbol,Symbol}, Array{Any,1}};
    stencil::Array{Array{Int,1},1}
    stencilBoundary::SparseVector{Int64,Int64}
    type::Symbol;
    kubWeights::Array{Float64,2};
    kubPoints::Array{Float64,2};
    taskRecovery::Bool;
    advection::Bool;
    recoveryOrders::Tuple;
end


#Konstruktoren

function femProblem(m::mesh, femType::Dict{Symbol, Array{Symbol,1}};stencilOrder=0,advection::Bool=true, taskRecovery::Bool=false, t::Symbol=:boussinesq, g::Int64=9, compoundMethod::Symbol=:none)
    sol=Dict{Float64, solution}()
    kubPoints, kubWeights=getKub(g, m.meshType);
    dF=Dict{Symbol, degF}()

    ordEdgesB, nebP, nebC=getOrderBoundary(m.boundaryEdges);
    ordVerticesB, nvbP, nvbC=getOrderBoundary(m.boundaryVertices);

    femElements=Set{Symbol}()
    recoveryType=Dict{Symbol, Int}()
    for k in keys(femType)
        if length(femType[k])==3
            recoveryType[k]=1
            nFEM=3
        else
            recoveryType[k]=0
            nFEM=length(femType[k])
        end
        for s in 1:nFEM
            push!(femElements,femType[k][s]);
        end
    end
    if t==:boussinesq
        recoveryOrders=(recoveryType[:p],recoveryType[:b],recoveryType[:v])
    elseif t==:compressible
        recoveryOrders=(recoveryType[:rhoTheta],recoveryType[:rhoV],recoveryType[:rho])
    elseif t==:shallow
        recoveryOrders=(recoveryType[:h],recoveryType[:hV])
    elseif t==:linshallow
        recoveryOrders=(recoveryType[:h],recoveryType[:v])
    end

    for k in femElements
        dF[k]=degF(m,k,ordEdgesB,nebP,nebC,ordVerticesB,nvbP,nvbC,kubPoints);
    end

    if compoundMethod!=:none
        compoundData=createCompoundData(compoundMethod,femElements,m);
    else
        compoundData=createCompoundData();
    end

    if taskRecovery && !iszero(stencilOrder)
        stencil, stencilBoundary=getStencil(m,stencilOrder)
    else
        stencil=Array{Array{Int,1},1}();
        stencilBoundary=spzeros(Int,0)
    end
    edgeData=Array{Array{Int64,1},1}();
    massM=Dict();
    massMB=Dict();
    stiffM=Dict();
    recoveryM=Dict();
    loadV=Dict();
    bV=Dict();
    s=Set{Symbol}([:poisson,:boussinesq,:compressible,:shallow,:linshallow]);
    !in(t,s) && error("Die Methode $t ist keine zulässige Eingabe. Möglich sind $s");
    femProblem(m,bV,dF,femType,edgeData,compoundData,sol,massM,massMB,stiffM,recoveryM,stencil,stencilBoundary,t,kubWeights, kubPoints, taskRecovery, advection, recoveryOrders);
end

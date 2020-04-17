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
    type::Symbol;
    kubWeights::Array{Float64,2};
    kubPoints::Array{Float64,2};
    taskRecovery::Bool;
    advection::Bool;
end


#Konstruktoren

function femProblem(m::mesh, femType::Dict{Symbol, Array{Symbol,1}};advection::Bool=true, taskRecovery::Bool=false, t::Symbol=:boussinesq, g::Int64=9, compoundMethod::Symbol=:none)
    sol=Dict{Float64, solution}()
    kubPoints, kubWeights=getKub(g, m.meshType);
    dF=Dict{Symbol, degF}()

    ordEdgesB, nebP, nebC=getOrderBoundary(m.boundaryEdges);
    ordVerticesB, nvbP, nvbC=getOrderBoundary(m.boundaryVertices);

    femElements=Set{Symbol}()
    for k in collect(keys(femType))
        for s in femType[k]
            push!(femElements,s);
        end
    end

    for k in femElements
        dF[k]=degF(m,k,ordEdgesB,nebP,nebC,ordVerticesB,nvbP,nvbC,kubPoints);
    end

    if compoundMethod!=:none
        compoundData=createCompoundData(compoundMethod,femElements,m);
    else
        compoundData=createCompoundData();
    end

    edgeData=Array{Array{Int64,1},1}();
    massM=Dict();
    massMB=Dict();
    stiffM=Dict();
    loadV=Dict();
    bV=Dict();
    s=Set{Symbol}([:poisson,:boussinesq,:compressible,:shallow]);
    !in(t,s) && error("Die Methode $t ist keine zulässige Eingabe. Möglich sind $s");
    femProblem(m,bV,dF,femType,edgeData,compoundData,sol,massM,massMB,stiffM,t,kubWeights, kubPoints, taskRecovery, advection);
end

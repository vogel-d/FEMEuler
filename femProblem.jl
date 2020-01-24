mutable struct femProblem
    mesh::mesh;
    boundaryValues::Dict{Tuple{Symbol,Symbol}, Array{AbstractFloat,1}};
    degFBoundary::Dict{Symbol, degF};
    femType::Dict{Symbol, Array{Symbol,1}};
    edgeData::Array{Array{Int,1},1};
    solution::Dict{AbstractFloat, solution};
    massM::Dict{Symbol, SuiteSparse.UMFPACK.UmfpackLU{Float64,Int}};
    massMBoundary::Dict{Symbol, SuiteSparse.UMFPACK.UmfpackLU{Float64,Int}};
    stiffM::Dict{Symbol, SparseMatrixCSC{AbstractFloat,Int}};
    type::Symbol;
    kubWeights::Array{AbstractFloat,2};
    kubPoints::Array{AbstractFloat,2};
    taskRecovery::Bool;
    advection::Bool;
end


#Konstruktoren

function femProblem(m::mesh, femType::Dict{Symbol, Array{Symbol,1}};advection::Bool=true, taskRecovery::Bool=false, t::Symbol=:boussinesq, g::Int=9)
    sol=Dict{AbstractFloat, solution}()
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
    edgeData=Array{Array{Int,1},1}();
    massM=Dict();
    massMB=Dict();
    stiffM=Dict();
    loadV=Dict();
    bV=Dict();
    s=Set{Symbol}([:poisson,:boussinesq,:compressible]);
    !in(t,s) && error("Die Methode $t ist keine zulässige Eingabe. Möglich sind $s");
    femProblem(m,bV,dF,femType,edgeData,sol,massM,massMB,stiffM,t,kubWeights, kubPoints, taskRecovery, advection);
end

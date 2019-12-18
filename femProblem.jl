mutable struct femProblem
    mesh::mesh;
    boundaryValues::Dict{Tuple{Symbol,Symbol}, Array{Float64,1}};
    degFBoundary::Dict{Symbol, degF};
    femType::Dict{Symbol, Array{Symbol,1}};
    edgeData::Array{Array{Int64,1},1};
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

function femProblem(m::mesh, femType::Dict{Symbol, Array{Symbol,1}};advection::Bool=true, taskRecovery::Bool=false, t::Symbol=:boussinesq, g::Int64=9)
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
    edgeData=Array{Array{Int64,1},1}();
    massM=Dict();
    massMB=Dict();
    stiffM=Dict();
    loadV=Dict();
    bV=Dict();
    s=Set{Symbol}([:poisson,:boussinesq,:combressible]);
    !in(t,s) && error("Die Methode $t ist keine zulässige Eingabe. Möglich sind $s");
    femProblem(m,bV,dF,femType,edgeData,sol,massM,massMB,stiffM,t,kubWeights, kubPoints, taskRecovery, advection);
end
#=
function femProblem(meth::Symbol, nx::Int64, ny::Int64, femType::Dict{Symbol,Array{Symbol,1}}, cond::Tuple{Symbol,Symbol};advection::Bool=true,  taskRecovery::Bool=false,  t::Symbol=:boussinesq,
                    g::Int64=9, xl::Float64=0.0, yl::Float64=0.0,xr::Float64=Float64(nx), yr::Float64=Float64(ny))
    if meth==:quad
        m=generateRectMesh(nx,ny,cond[1],cond[2],xl,xr,yl,yr);
    elseif meth==:qtri
        m=generateTriMesh(nx,ny,xl,yl,xr,yr);
    elseif meth==:tri
        m=generateTriMesh2(nx,ny,xl,yl,xr,yr);
    elseif meth==:hex
        error("Hexagonale Gitter können nur generiert und noch nicht gelöst werden.")
    else
        error("Keine zulässige Gitterstruktur! Mögliche Eingaben: {:quad, :qtri, :tri}")
    end

    kubPoints, kubWeights=getKub(g, m.meshType);
    sol=Dict{Float64, solution}();
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
    #=
    femElements=collect(femElements)
    h = @distributed (append!) for k in 1:length(cfemElements)
        [degF(m,cfemElements[k],boundary,b,kubPoints)];
    end
    for k in 1:length(cfemElements)
        dF[cfemElements[k]]=h[k]
    end
    =#
    edgeData=Array{Array{Int64,1},1}();
    massM=Dict();
    massMB=Dict();
    stiffM=Dict();
    loadV=Dict();
    bV=Dict();
    s=Set{Symbol}([:poisson,:boussinesq,:compressible]);
    !in(t,s) && error("Die Methode $t ist keine zulässige Eingabe. Möglich sind $s");
    femProblem(m,bV,dF,femType,edgeData,sol,massM,massMB,stiffM,t,kubWeights, kubPoints, taskRecovery, advection);
end
=#

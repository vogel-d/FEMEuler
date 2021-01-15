using HDF5, JLD
function saveProblem(p::femProblem, filename::String)
    save((@__DIR__)*"/JLD/"*filename*".jld", "femProblem", p)
    return nothing;
end
function loadProblem(filename::String)
    return load((@__DIR__)*"/JLD/"*filename*".jld")["femProblem"];
end

function saveMesh(m::mesh, filename::String)
    save((@__DIR__)*"/JLD/"*filename*".jld", "mesh", m)
    return nothing;
end
function loadMesh(filename::String)
    return load((@__DIR__)*"/JLD/"*filename*".jld")["mesh"];
end

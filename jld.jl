using HDF5, JLD
function saveProblem(p::femProblem, filename::String)
    save((@__DIR__)*"/VTK/"*filename*".jld", "femProblem", p)
    return nothing;
end
function loadProblem(filename::String)
    return load((@__DIR__)*"/VTK/"*filename*".jld")["femProblem"];
end

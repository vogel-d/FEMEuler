function getEdgeType(n::Array{AbstractFloat,1})
    if n==[0.0,-1.0]
        edgeType=1;
    elseif n==[-1.0,0.0]
        edgeType=2;
    elseif n==[0.0, 1.0]
        edgeType=3;
    elseif n==[1.0,0.0]
        edgeType=4;
    elseif n==[0.7071067811865475244, 0.7071067811865475244]
        edgeType=1;
    else
        error("Kein zulässiger Normalenvektor!");
    end

    return edgeType
end

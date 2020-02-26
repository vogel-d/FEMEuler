function getEdgeType(n::Array{Float64,1})
    if n==[0.0,-1.0]
        edgeType=1;
    elseif n==[-1.0,0.0]
        edgeType=2;
    elseif n==[0.0, 1.0]
        edgeType=3;
    elseif n==[1.0,0.0]
        edgeType=4;
    elseif isequal(n,[1/sqrt(2), 1/sqrt(2)])
        edgeType=3;
    else
        error("Kein zulässiger Normalenvektor: $(n)");
    end

    return edgeType
end
